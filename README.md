# ParamSim

The parametrized simulation is currently a module in GArSoft v02_04_00 and is complementary to GArSoft's full reconstruction approach. GArSoft is Multipurpose Detector's (MPD) official simulation and analysis software framework. v02_04_00 of GArSoft is the version tagged for CDR studies. 

Two reconstruction approaches are taken to parametrize the detector response: 

1) Gluckstern relation (https://www.sciencedirect.com/science/article/pii/S0168900209009474) for reconstructing long tracks curving up in magnetic field  
2) Particle range for reconstructing short tracks that stop in the High pressure gas argon TPC component of the MPD and do not curve up in magnetic field

In addition to reconstructing tracks, a dE/dx-based PID is implemented in the module. This is based on Tom Junk's parametrization of PEP-4's dE/dx curve: https://home.fnal.gov/~trj/mpd/dedx_sep2019/

The module is designed to take GArSoft's analysis tree, anatree as input and produce a so called "cafanatree" ntuple as output. A description of cafanatree tree variables are as the following: 

* Event: an art event (not neutrino-interaction event)

* Run: run #

* SubRun: subrun #

- Generator-level information:

  * mode: running mode (neutrino vs. anti-neutrino)

  * q2: momentum transfer

  * w: hadronic invariant mass

  * y: inelasticity

  * x: Bjorken-x variable

  * theta: angle between incoming and outgoing lepton

  * t: a variable used in coherent pion analysis

  * ntype: neutrino type/flavor 

  * ccnc: charged current or neutral current

  * gint: a variable used in coherent pion analysis

  * tgtpdg: target PDG (nucleus-level variable)

  * weight: number of protons-on-target

  * gt_t: a variable used in coherent pion analysis

  * intert: interaction type 

  * mcnupx: x-component momentum of neutrino

  * mcnupy: y-component momentum of neutrino

  * mcnupz: z-component momentum of neutrino

  * vertx: vertex of primary mc particle in x

  * verty: vertex of primary mc particle in y

  * vertz: vertex of primary mc particle in z
