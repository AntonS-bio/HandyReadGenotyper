import pickle

class Model:
    def __init__(self, name: str) -> None:
        self.name=name
        self._sensitivity=-1.0
        self._specificity=-1.0
        self._training_samples=-1
        self._testing_samples=-1
        self._testing_samples_misclassfied=-1
        self._nucletoide_seq=""

    @property
    def sensitivity(self) -> float:
        return self._sensitivity

    @sensitivity.setter
    def sensitivity(self, value: float):
        self._sensitivity = value

    @property
    def specificity(self) -> float:
        return self._specificity

    @specificity.setter
    def specificity(self, value: float):
        self._specificity = value

    @property
    def training_samples(self) -> int:
        return self._training_samples

    @training_samples.setter
    def training_samples(self, value: int):
        self._training_samples = value

    @property
    def testing_samples(self) -> int:
        return self._testing_samples

    @testing_samples.setter
    def testing_samples(self, value: int):
        self._testing_samples = value

    @property
    def testing_samples_misclassfied(self) -> int:
        return self._testing_samples_misclassfied

    @testing_samples_misclassfied.setter
    def testing_samples_misclassfied(self, value: int):
        self._testing_samples_misclassfied = value

    @property
    def nucletoide_seq(self) -> float:
        return self._nucletoide_seq

    @nucletoide_seq.setter
    def nucletoide_seq(self, value: float):
        self._nucletoide_seq = value
