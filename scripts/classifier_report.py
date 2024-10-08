import pickle
from typing import List, Dict
from read_classifier import ClassificationResult, Classifier, ModelsData
from os.path import expanduser, basename
from pathlib import Path
from datetime import datetime
import json
from collections import Counter
from map import MappingResult

OPENER='<tr><td>'
OPENER_BOLD='<tr class=table_header><td>'
OPENER_ALT='<tr class=alt_row><td>'
CLOSER='</td></tr>'
MIDDLE='</td><td>'
MIDDLE_WIDE='</td class=td_alt_wide><td>'
POSITIVE_CASES_THRESHOLD=10
LOW_READ_WARNING=10

class ClasssifierReport:

    def __init__(self, output_file: str, models_file: str, positive_amplicons_threshold: float, sample_labels:Dict[str, str]={}, 
                 mapping_results: List[MappingResult]=[]) -> None:
        self._positive_amplicons_threshold = positive_amplicons_threshold/100

        self._model_file=models_file
        self.output_file=open( expanduser(output_file), "w" )
        self.genotypes={}
        with open(expanduser(models_file), "rb") as model_file:
            self.models_data: ModelsData = pickle.load(model_file)
            self.trained_models: Dict[str, Classifier]= self.models_data.classifiers #Dict[key=amplicon name, value=corresponding model]
        if "hierarchy" in self.models_data.metadata:
            self.genotypes = self.models_data.metadata["hierarchy"]
        self.sample_labels=sample_labels
        self.sample_hyperlinks={}
        self.mapping_results={}
        for item in mapping_results:
            self.mapping_results[item.sample_name] = item.total_reads

    def _clean_sample(self, sample_name) -> str:
        return Path(sample_name).stem

    def _write_main_header(self, file_handle):
        file_handle.write("<table>\n")
        file_handle.write("\t<tbody>\n")
        #header_values=["Sample","Has Target Organism*", "Implied Genotypes","Implied AMR", "Amplicons by # of unknown SNPs", "Sequencing Reads","Mapped, correct length", "Mapped, but too short/long"]#," Highest read count","Lowest read count"]
        header_values=["Sample", "Amplicons by # of unknown SNPs", "Implied Genotypes","Implied AMR", "Sequencing Reads","Mapped, correct length", "Mapped, but too short/long"]#," Highest read count","Lowest read count"]
        header_line=OPENER_BOLD+MIDDLE.join(header_values)+CLOSER
        file_handle.write("\t"+header_line)
        file_handle.write("\n")

    def _write_styles(self):
        with open( Path(__file__).with_name('html_head.txt') ,"r" ) as styles_data:
            for line in styles_data:
                self.output_file.write(line)
            self.output_file.write("\n")

    def _write_navigation(self, file_handle, locations: List[str]):
        file_handle.write("\t<nav>\n")
        for location in locations:
            file_handle.write('\t\t<a href="#top">'+self._clean_sample(location)+'</a>\n')
        file_handle.write("\t</nav>\n")

    def _get_sample_results(self, sample: str, results: List[ClassificationResult]) -> List[ClassificationResult]:
        return [f for f in results if f.sample==sample]

    def _insert_paragraph(self, number):
        for i in range(0,number):
            self.output_file.write("<br/>")
        self.output_file.write("\n")

    def _write_gt(self, key: str, values:Dict[str, Dict], depth:int):
        if depth!=0:
            self.genotype_row_labels[key]="<p>"+"".join(["&emsp;"]*(depth-1))+"&#8627;"+key+"</p>"
        else:
            self.genotype_row_labels[key]="<p>"+key+"</p>"
        self.genotype_row_order[key]=len(self.genotype_row_order)
        if "children" in values and len(values["children"])>0:
            for key, item in values["children"].items():
                self._write_gt( key, item, depth+1)

    def sample_display_name(self, sample:str):
        sample_cell_value=self._clean_sample(sample)
        if self._clean_sample(sample) in self.sample_labels:
            sample_cell_value=sample_cell_value+"<br>"+self.sample_labels[self._clean_sample(sample)]
        return sample_cell_value

    def get_most_recent_model_element(self):
        models_list=sorted(set([(f.model_timestamp) for f in self.results]))
        return models_list[-1].strftime("%d/%b/%y")

    def _unamplified_amr_target(self, result: ClassificationResult) -> str:
        #Get strikethrough names of the AMR amplicons that have not been amplified.
        is_amr=False
        if result.amplicon.name in self.trained_models:
            for snp in self.trained_models[result.amplicon.name].genotype_snps:
                if snp.is_amr:
                    is_amr = True
        if not is_amr:
            return ""
        if result.positive_cases == 0:
            return "<del>"+result.amplicon.name+"</del>"
        else:
            return ""

    def get_transient_amplicons(self, results: List[ClassificationResult]) -> str:
        values=[]
        for result in results:
            amplicon_name=result.amplicon.name
            if result.positive_cases>POSITIVE_CASES_THRESHOLD and self.trained_models[amplicon_name].amplicon_transient:
                values.append(amplicon_name)
        return values

    ####START Write main summary table
    def _write_summary(self):
        self.output_file.write('<div class=header_line><a name="Summary">Summary of results</div>\n')
        self.output_file.write('<div>Model name: '+basename(self._model_file).replace(".pkl","")+
                                ', model date: '+self.get_most_recent_model_element()+'</div>\n')
        self._insert_paragraph(1)

        self._write_main_header(self.output_file)

        self.genotype_row_labels={}
        self.genotype_row_order={}
        for key, item in self.genotypes.items():
            self._write_gt(key, item, 0)

        self.sample_gts_w_frequency: Dict[str, List[str]]={}
        self.sample_gts: Dict[str, List[str]]={}
        self.empty_samples=set()
        add_read_length_flag=False
        for i, sample in enumerate(self.samples):
            sample_results=self._get_sample_results(sample,self.results)
            mapped_reads=[f.positive_cases+f.negative_cases for f in sample_results]
            wrong_len_reads=sum([f.wrong_len_reads[f.amplicon.name] for f in sample_results])
            total_mapped_reads=sum(mapped_reads)
            if total_mapped_reads>0 and wrong_len_reads/total_mapped_reads > 3:
                total_mapped_reads=str(total_mapped_reads)+"^"
                add_read_length_flag=True
            highest_cov = sample_results[ mapped_reads.index(max(mapped_reads)) ]
            highest_cov=f'{len(highest_cov.predicted_classes)}<br>({highest_cov.amplicon.name})'
            lowest_cov=sample_results[ mapped_reads.index(min(mapped_reads)) ]
            lowest_cov=f'{len(lowest_cov.predicted_classes)}<br>({lowest_cov.amplicon.name})'
            #threshold settings: positive_amplicons_threshold% of amplicons must have >POSITIVE_CASES_THRESHOLD target reads
            has_target_org="Yes"
            amr_snps=""
            gts=""
            sample_snp_diagram=self.snp_frequency_diagram(sample)
            #check if the number of non-transient amplicons with sufficient number of reads is above the presence threshold
            non_transient_amplicons=sum([1 for f in sample_results if not self.trained_models[f.amplicon.name].amplicon_transient])
            non_transient_with_reads=sum( [1 for f in sample_results if f.positive_cases>POSITIVE_CASES_THRESHOLD and not self.trained_models[f.amplicon.name].amplicon_transient] )
            # if non_transient_with_reads>non_transient_amplicons*self._positive_amplicons_threshold:
            #     has_target_org="Yes"
            # else:
            #     self.empty_samples.add(self._clean_sample(sample))
            
            if has_target_org=="No": #do not show implied GTs if no target organism to avoid confusion.
                gts=""
            else:
                for target, option in zip( [self.sample_gts, self.sample_gts_w_frequency], [False, True]  ):
                    gts=sorted(set([f.calculate_genotype(self.trained_models[f.amplicon.name].genotype_snps, LOW_READ_WARNING, option) for f in sample_results]))
                    #calculate_genotype( ) ouputs Tuple[str, bool] with SNP name (gt) and bool to indicate AMR (True) or not (False)
                    if "" in gts:
                        gts.pop(gts.index(""))
                    target[ self._clean_sample(sample) ]=set([f[0] for f in gts])

                missing_amr_products=[f for f in [self._unamplified_amr_target(result) for result in sample_results] if f!=""]
                amr_snps=", ".join([f[0] for f in  gts if f[1] and f[0]!=""] + missing_amr_products + self.get_transient_amplicons(sample_results) )
                gts=", ".join([f[0] for f in  gts if not f[1] and f[0]!=""])
            sample_cell_value=self.sample_display_name(sample)
            if self._clean_sample(sample) in self.mapping_results: #only works if mapping was done by HRG
                sequenced_reads=self.mapping_results[self._clean_sample(sample)]
            else:
                sequenced_reads="Unavailable"
            sample_line_values=['<a href="#'+self._clean_sample(sample)+'">'+sample_cell_value+"</a>",
                                #has_target_org,
                                sample_snp_diagram,
                                gts,
                                amr_snps,
                                #"<td class=diagram>"+sample_snp_diagram,
                                sequenced_reads,
                                total_mapped_reads,
                                wrong_len_reads]
                                #highest_cov,
                                #lowest_cov ]
            if i % 2 == 0: #Alternate row colours
                table_line=OPENER_ALT+MIDDLE.join([str(f) for f in sample_line_values])+CLOSER
            else:
                table_line=OPENER+MIDDLE.join([str(f) for f in sample_line_values])+CLOSER
            #this is a hack to accomodate special formatting for snp frequency diagram cell
            table_line=table_line.replace("<td><td class=diagram>","<td class=diagram>")
            self.output_file.write("\t\t"+table_line+"\n")

        self.output_file.write("\t</tbody>\n")
        self.output_file.write("</table>\n")
        #self.output_file.write(f'<p>* Defined as at least {self._positive_amplicons_threshold*100}% of amplicons that <a href="#ModelSignatures">are part of core genome</a> having at least {POSITIVE_CASES_THRESHOLD} target organism reads</p>')
        self.output_file.write(f'<p>(?) Prediction is based on amplicon with fewer than {LOW_READ_WARNING} reads</p>')
        self.output_file.write(f'<p><del>AMR Amplicon</del> means no reads were found for this AMR amplicon</p>')
        if add_read_length_flag:
            self.output_file.write(f'<p>^ The number of reads either longer or shorter than amplicon is >3x higher than reads that match amplicon length. If this is unexpected, consider using options "-s" or "-l". </p>')
        self.output_file.write('<p>Report generated on '+datetime.now().strftime("%d/%b/%y at %H:%M")+'. <a href="#ModelSignatures">See models list</a></p>')

        self._insert_paragraph( 2)
    ####END Write main summary table

    ####START Write models' summary
    def _write_model_summary(self):
        self.output_file.write('<div class=header_line><a name="ModelSignatures">Model Signatures</a><div><a href="#Summary">Back to Summary</a></div></div>\n')
        self._insert_paragraph(1)
        self.output_file.write("<table>\n")
        self.output_file.write("\t<tbody>\n")
        header_values=["Amplicon", "Is Core","Signature", "Date Trained"]
        table_header=OPENER_BOLD+MIDDLE.join(header_values)+CLOSER
        self.output_file.write(table_header+"\n")
        models_list=sorted(set([(f.amplicon.name, f.model_fingerprint,f.model_timestamp) for f in self.results]))
        for i, result in enumerate(models_list):
            is_core=not self.trained_models[result[0]].amplicon_transient
            row_values=[result[0], is_core, result[1], result[2].strftime("%d/%b/%y at %H:%M")]
            if i % 2 == 0: #Alternate row colours
                table_row=OPENER_ALT+MIDDLE.join([str(f) for f in row_values])+CLOSER
            else:
                table_row=OPENER+MIDDLE.join([str(f) for f in row_values])+CLOSER
            self.output_file.write("\t\t"+table_row+"\n")
        self.output_file.write("<table>\n")
        self.output_file.write("\t<tbody>\n")
        self.output_file.write("</div>\n")
    ####END Write models' summary

    ####START Write classification results
    def _write_classification_results(self):
        amplicon_names=sorted(list(set([f.amplicon.name for f in self.results])))
        self.output_file.write("<div class=sample_reports>\n")
        for sample in self.samples:
            sample_cell_value=self.sample_display_name(sample)
            self.output_file.write('<div class=header_line><a name="'+self._clean_sample(sample)+'">'+sample_cell_value+'</a><div><a href="#Summary">Back to Summary</a></div></div>\n')
            self._insert_paragraph(1)
            self.output_file.write("<table>\n")
            self.output_file.write("\t<tbody>\n")

            header_values=["Amplicon","Tot. Reads", "% Target Org.", "# Target Org.","SNPs", "Implied Genotypes", "Implied AMR","Unknown SNPs"]
            table_header=OPENER_BOLD+MIDDLE.join(header_values)+CLOSER
            self.output_file.write("\t"+table_header+"\n")
            sample_results=self._get_sample_results(sample,self.results)
            for i, amplicon_name in enumerate(amplicon_names):
                untrained_model=self.trained_models[amplicon_name]._not_trained
                untrained_suffix="^" if untrained_model else ""
                amplicon_results=[f for f in sample_results if f.amplicon.name==amplicon_name]
                if len(amplicon_results)==0:
                    raise ValueError(f'Amplicon {amplicon_name} is missing for sample {sample}')
                amplicon_result: ClassificationResult = amplicon_results[0]
                first_cell='<a href="#'+self._clean_sample(sample)+'_'+amplicon_name+'_consensus"/>'+amplicon_name
                if len(amplicon_result.predicted_classes)<=POSITIVE_CASES_THRESHOLD:
                    row_values=[first_cell,"*", "-","-","-", "-","-", "-"]
                else:
                    tot_reads=len(amplicon_result.predicted_classes)
                    target_org_reads_share='{0:.1%}'.format(amplicon_result.positive_cases/tot_reads)+untrained_suffix
                    target_org_reads_count=str(amplicon_result.positive_cases)+untrained_suffix
                    snp_num=amplicon_result.num_mismatches
                    unknown_snps=amplicon_result.num_unknown_snps(self.trained_models[amplicon_name].genotype_snps)
                    known_snps=[f for f in amplicon_result.calculate_genotype(self.trained_models[amplicon_name].genotype_snps, LOW_READ_WARNING)]
                    known_gts=known_snps[0] if known_snps[1] is False else ""
                    known_amrs=known_snps[0] if known_snps[1] is True else ""
                    row_values=[first_cell + untrained_suffix, tot_reads,target_org_reads_share,target_org_reads_count,snp_num, known_gts, known_amrs,unknown_snps]
                if i % 2 == 0: #Alternate row colours
                    table_row=OPENER_ALT+MIDDLE.join([str(f) for f in row_values])+CLOSER
                else:
                    table_row=OPENER+MIDDLE.join([str(f) for f in row_values])+CLOSER
                self.output_file.write("\t\t"+table_row+"\n")
            self.output_file.write("\t<tbody>\n")
            self.output_file.write("<table>\n")
            self.output_file.write(f'* Amplicon has fewer than {POSITIVE_CASES_THRESHOLD+1} total reads<br>')
            self.output_file.write(f'^ Using model that was trained without negative cases - results are likely to be less reliable')
            self._insert_paragraph(5)
    ####END Write classification results

    ####START write genotypes support tables
    def _write_gt_suppport(self):
        if len(self.genotypes)==0:
            #either genotype file not provided, or it's empty.
            self.output_file.write("<div class=header_line>Genotypes summary:<br>model file does not contain genotypes hierarchy</div>\n")
            self._insert_paragraph(2)
            return

        self.output_file.write("<div class=header_line>Genotypes summary</div>\n")
        self._insert_paragraph(1)
        self.output_file.write("<table>\n")
        self.output_file.write("\t<tbody>\n")
        header_values=["Genotype"]+[self.sample_display_name(f) for f in self.samples]
        table_header=OPENER_BOLD+MIDDLE.join(header_values)+CLOSER+"\n"
        self.output_file.write(table_header)
        sorted_genotypes_labels=sorted(self.genotype_row_labels, key= lambda x: self.genotype_row_order[x])
        local_opener='<td class="td_alt">'
        for i, label in enumerate(sorted_genotypes_labels):
            row_values=[ self.genotype_row_labels[label] ]
            for sample in [self._clean_sample(f) for f in self.samples]:
                # if sample in self.empty_samples:
                #     row_values.append("n.d.")
                # else:
                row_values.append( "S" if  label in self.sample_gts[sample] else "" )
            if i % 2 == 0: #Alternate row colours
                table_row=OPENER_ALT.replace("<td>",local_opener) + MIDDLE.join([str(f) for f in row_values])+CLOSER
            else:
                table_row=OPENER.replace("<td>",local_opener) + MIDDLE.join([str(f) for f in row_values])+CLOSER
            self.output_file.write("\t\t"+table_row+"\n")

        self.output_file.write("\t</tbody>\n")
        self.output_file.write("</table>\n")
        self.output_file.write("<p>S - reads support this genotype;<br>n.d. - target organism not detected in the sample.")
        self.output_file.write("</div>\n")
        self._insert_paragraph( 1)
    ####END write genotypes support tables

    def snp_frequency_diagram(self, sample: str):
        HTML_SPACE = "&nbsp;"
        amplicon_names=sorted(list(set([f.amplicon.name for f in self.results])))
        sample_results=self._get_sample_results(sample,self.results)
        unknown_snps=[]
        for i, amplicon_name in enumerate(amplicon_names):
            amplicon_results=[f for f in sample_results if f.amplicon.name==amplicon_name]
            if len(amplicon_results)==0:
                raise ValueError(f'Amplicon {amplicon_name} is missing for sample {sample}')
            amplicon_result: ClassificationResult = amplicon_results[0]
            if amplicon_result.positive_cases>POSITIVE_CASES_THRESHOLD:
                unknown_snps.append(amplicon_result.num_unknown_snps(self.trained_models[amplicon_name].genotype_snps))
        unknown_snps=Counter(unknown_snps)
        diagram=""
        for x in range(0, 5):
            if x in unknown_snps:
                diagram=diagram+str(unknown_snps[x])+","+HTML_SPACE
            else:
                diagram=diagram+"0,"+HTML_SPACE
        over_five=sum([value for key, value in unknown_snps.items() if key>=5])
        diagram=diagram+str(over_five)
        return diagram
        # diagram="<br/><br/>"+"0"+HTML_SPACE+HTML_SPACE+"3"+HTML_SPACE+">5"
        # for y in range(1, 6):
        #     new_line=""
        #     for x in range(0, 6):
        #         new_line+="-" if x in unknown_snps and unknown_snps[x]>=y else HTML_SPACE
        #     over_fives=sum ([value for key, value in unknown_snps.items() if key>5] )
        #     new_line+="-" if over_fives>=y else HTML_SPACE
        #     diagram=new_line+"<br/>"+diagram
        return diagram

    def _write_consensus_sequences(self):
        self._insert_paragraph(1)
        self.output_file.write("<div class=header_line>Consensus sequences</div>\n")
        self._insert_paragraph(1)
        self.output_file.write("<table>\n")
        self.output_file.write("\t<tbody>\n")
        header_values=["Sample", "Amplicon", "# of target reads", "Mismatches <br/> vs Reference", "Consensus sequence"]
        table_header=OPENER_BOLD+MIDDLE.join(header_values)+CLOSER+"\n"
        self.output_file.write(table_header)

        for i, result in enumerate(self.results):
            snp_values: List[str] = []
            for known_snp in result.known_snps(self.trained_models[result.amplicon.name].genotype_snps):
                allele_frequency=" ("+'{:.0%}'.format(result.consensus_frequency[known_snp["Pos"]])+") "
                snp_value="Pos:"+str(known_snp["Pos"]) + allele_frequency
                snp_value+=", Change:"+result.amplicon.ref_seq.sequence[known_snp["Pos"]] + ">"+known_snp["Nucleotide"]
                snp_value+=", Implication: "+known_snp["SNP"].get_genotype(known_snp["Nucleotide"])
                snp_values.append(snp_value)
            for known_snp in result.unknown_snps(self.trained_models[result.amplicon.name].genotype_snps):
                allele_frequency=" ("+'{:.0%}'.format(result.consensus_frequency[known_snp["Pos"]])+") "
                snp_value="Pos:"+str(known_snp["Pos"]) + allele_frequency
                snp_value+=", Change:"+result.amplicon.ref_seq.sequence[known_snp["Pos"]] + ">"+known_snp["Nucleotide"]
                snp_value+=", Implication: Unknown"
                snp_values.append(snp_value)
            snp_values='<br/>'.join(snp_values)
            short_name=self._clean_sample(result.sample)
            hyperlink_name=short_name+"_"+result.amplicon.name+"_consensus"
            first_cell='<a name="'+hyperlink_name+'" href="#'+short_name+'">'+self.sample_display_name(result.sample)+'</a>'
            row_values=[first_cell, result.amplicon.name, result.positive_cases, snp_values, result.consensus]
            if i % 2 == 0: #Alternate row colours
                table_row=OPENER_ALT + MIDDLE_WIDE.join([str(f) for f in row_values])+CLOSER
            else:
                table_row=OPENER + MIDDLE_WIDE.join([str(f) for f in row_values])+CLOSER
            self.output_file.write("\t\t"+table_row+"\n")



        self.output_file.write("\t</tbody>\n")
        self.output_file.write("</table>\n")
        self.output_file.write("</div>\n")

    def create_report(self, results: List[ClassificationResult]):
        self.results: List[ClassificationResult] = results
        self.samples: List[str]=sorted(set([f.sample for f in self.results]))

        self._write_styles()
        self.output_file.write("<html>\n")
        self.output_file.write("<body>\n")

        self._write_summary()

        self._write_gt_suppport()

        self._write_classification_results()

        self._write_model_summary()

        self._write_consensus_sequences()
        self.output_file.write("</body>\n")
        self.output_file.write("</html>\n")
        self.output_file.close()
