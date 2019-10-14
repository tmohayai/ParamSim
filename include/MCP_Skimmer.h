#ifndef MCP_SKIMMER_H
#define MCP_SKIMMER_H 1

#include <string>

#include "TFile.h"
#include "TTree.h"

class MCP_Skimmer
{
public:

    /* Default Constructor */
    MCP_Skimmer();

    /* Constructor */
    MCP_Skimmer(std::string infilename, std::string outfilename);

    /* Default Copy Constructor */
    MCP_Skimmer(const MCP_Skimmer &) = default;

    /* Destructor */
    ~MCP_Skimmer();

    /* Method to skim the MCPs from the TFile */
    void SkimMCParticle();

    /* Method to book the output TFile */
    bool BookTFile();

    /* Method to fill the TTree */
    void FillTTree();

    /* Method to write the TTree */
    void WriteTTree();

    /* Method to close the TFile */
    void CloseTFile();

    /* Method to clear the event variables */
    void ClearVectors();

protected:

private:
    std::string _outfile;
    TFile* _skimfile;
    TTree* _skimtree;

    std::string _infile;
    TFile* _intfile;
    TTree* _inttree;

    int           	 _Event;
    int           	 _SubRun;
    int           	 _Run;

    std::vector<float> _MC_Q2;
    std::vector<float> _MC_W;
    std::vector<float> _MC_Y;
    std::vector<float> _MC_X;
    std::vector<float> _MC_Theta;
    std::vector<float> _MC_T;
    std::vector<float> _MCVertX;
    std::vector<float> _MCVertY;
    std::vector<float> _MCVertZ;
    std::vector<float> _MCNuPx;
    std::vector<float> _MCNuPy;
    std::vector<float> _MCNuPz;
    std::vector<int>   _NType;
    std::vector<int>   _CCNC;
    std::vector<int>   _Mode;
    std::vector<int>   _Gint;
    std::vector<int>   _TgtPDG;
    std::vector<int>   _GT_T;
    std::vector<int>   _InterT;
    std::vector<float> _Weight;

    std::vector<int>     _PDG;
    std::vector<int>     _Mother;
    std::vector<int>     _PDGMother;
    std::vector<int>     _MCPTrkID;
    std::vector<float>   _MCPStartX;
    std::vector<float>   _MCPStartY;
    std::vector<float>   _MCPStartZ;
    std::vector<float>   _MCPStartPX;
    std::vector<float>   _MCPStartPY;
    std::vector<float>   _MCPStartPZ;
    std::vector<std::string>   _MCPProc;
    std::vector<std::string>   _MCPEndProc;
    std::vector<float>   _MCPEndX;
    std::vector<float>   _MCPEndY;
    std::vector<float>   _MCPEndZ;
    std::vector<float>   _TrajMCPX;
    std::vector<float>   _TrajMCPY;
    std::vector<float>   _TrajMCPZ;
    std::vector<int>     _TrajMCPTrajIndex;
};

#endif /* MCP_SKIMMER_H */
