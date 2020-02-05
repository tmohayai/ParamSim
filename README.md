# ParamSim

The parametrized simulation is currently a module in GArSoft v02_04_00 and is complementary to GArSoft's full reconstruction approach. GArSoft is Multipurpose Detector's (MPD) official simulation and analysis software framework (https://cdcvs.fnal.gov/redmine/projects/garsoft/wiki), and v02_04_00 of GArSoft is the version tagged for CDR studies. 

Following is the reconstruction approach taken to parametrize the detector response: 

1) Gluckstern relation (https://www.sciencedirect.com/science/article/pii/S0168900209009474) for reconstructing long tracks curling up in magnetic field  
2) Particle range for reconstructing short tracks that stop in the high pressure gas argon TPC (HPgTPC) component of the MPD and do not curking up in magnetic field

In addition to reconstructing tracks, a dE/dx-based PID is implemented in the module. This is based on Tom Junk's parametrization of PEP-4's dE/dx curve: https://home.fnal.gov/~trj/mpd/dedx_sep2019/ (PID matrices are in pid.root file)

The module is designed to take GArSoft's analysis tree, anatree as input and produce a so called "cafanatree" ntuple as output. A description of cafanatree tree variables are as the following: 

* Event: an art event (not neutrino-interaction event)

* Run: run #

* SubRun: subrun #

- Generator-level information (more information about generator-level variables can be found at: https://internal.dunescience.org/doxygen/classsimb_1_1MCNeutrino.html)

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

  * mcnupx: x-component momentum of neutrino in GeV/c

  * mcnupy: y-component momentum of neutrino in GeV/c

  * mcnupz: z-component momentum of neutrino in GeV/c

  * vertx: vertex of primary mc particle in x in cm

  * verty: vertex of primary mc particle in y in cm 

  * vertz: vertex of primary mc particle in z in cm
  
- GEANT4-level information: 

  * nFSP: number of mc particles -- note that this variable is not limited to "primary" mc particles that emerge from the neutrino interaction vertices, rather it includes all mc particles
 
  * pdgmother: pdg code of the particle that created the particle under consideration
 
  * mctrkid: track number of the particle that created the track under study

  * mctime: detector response time/time information of the particle with respect to neutrino interaction time (which is 0 nano seconds)

  * MCPStartX: starting position of the particle along the x-direction (starting position of the track) in cm 

  * MCPStartY: starting position of the particle along the y-direction (starting position of the track) in cm 

  * MCPStartZ: starting position of the particle along the z-direction (starting position of the track) in cm

  * truepx: truth-level x-component momentum of the particle in GeV/c
 
  * truepy: truth-level y-component momentum of the particle in GeV/c

  * truepz: truth-level z-component momentum of the particle in GeV/c

  * MCPEndX: end position of the particle along the x-direction (end position of the track) in cm 

  * MCPEndY: end position of the particle along the y-direction (end position of the track) in cm 

  * MCPEndZ: end position of the particle along the z-direction (end position of the track) in cm 

  * MCProc: a vector containing a string of the GEANT4 process that created a particle (e.g. for particles emerging from a neutrino interacion, the GEANT4 process is "primary")

  * MCEndProc: a vector containing a string of the GEANT4 end process of a particle (e.g. for e+e- pairs emerging from a photon conversion, the end process is "conv")

  * angle: truth-level angle with respect to the beam direction (along the z) 
 
  * truep: truth-level momentum of the particle in GeV/c

  * truepdg: truth-level PDG code of the particle 
  
- Reconstruction-level information (note that this is different from GArSoft reconstruction -- see above for more information about the parametrized reconstruction in the ParamSim module): 

  * recopid: reconstructed PDG code of the particle using the parametrized dE/dx

  * trkLen: length of the track

  * trkLenPerp: length of the track in the direction transverse to the beam 

  * preco: reconstructed particle momentum 

  * anglereco: reconstructed angle of the particle with respect to the beam direction (wrt the z-direction) 

  * erecon: reconstructed energy (so far only contains neutral particle energy reconstruction using ECAL) 

  * prob_arr: array of PID scores/probabilities (using log-likelihood) for reconstructed particles 
  
Check out the test directory for an example macro on how to read the cafanatree analysis ntuples that are produced as outputs of running the ParamSim module.   
