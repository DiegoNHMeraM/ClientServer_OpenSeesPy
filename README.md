# Client-Server application for vRTHS

## Description

The example provided here is a three-story planar frame with two bays is chosen as reference structure, with a total of 30 degrees of freedom (9 lateral, 9 vertical, and 12 rotational). The reference structure is divided into the numerical substructure (NS) and experimental substructure (ES). One column of the first story is taken as the ES, while the rest is modeled as NS. Since it is intended to work only with a uniaxial actuator, the boundary degrees of freedom are set to be only horizontal displacements.

It is important to mention that these substructures, numerical and experimental, can be replaced by the substructures desired by the user.

<img src="Figures/Subestructures.png" alt="Reference Structure" width="800"/>