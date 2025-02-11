from dataclasses import dataclass


@dataclass
class IonClass:
    name: str = ""
    Cout: float = 1.0
    Cin: float = 1.0
    valence: float = 1.0
    Erev: float = 0.0


# we'd love to be able to add ions later.. neuronsim doesn't support this yet
all_ions = [
    IonClass(name='Na', Cout=140.0, Cin=8.0, valence=+1),
    IonClass(name='K', Cout=4., Cin=140., valence=+1),
    IonClass(name='Cl', Cout=140., Cin=20., valence=-1),
]
