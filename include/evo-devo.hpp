#pragma once
#include <cstdint>
#include <unordered_map>

const int fieldsNumber      = 16;   
const int stemCellsTypes    = 32;

typedef int32_t RelativeStemCellIndex; 

enum Direction : uint8_t{
    RIGHT       = 0,    //0000
    LEFT        = 1,    //0001
    UP          = 2,    //0010
    DOWN        = 3,    //0011
    FORWARD     = 4,    //0100
    BACK        = 5,    //0101
    OUTWARDS    = 8,    //1000
    INWARDS     = 9     //1001
};


//The genome completely defines the behaviour of all stemCells during development
struct Genome_t{
    //Each autosome contains the behaviour of a different stemCells type
    struct {
        //Each gene contains the amplitude that the stemCell generates for its field, which type it will transition to if exposed to this field, and the threshold above/below which the transition happens, the direction to spawn a new stemCell, and the threshold associated with it
        struct {
            int8_t      pulseThreshold;     // threshold to generate a field pulse
            uint8_t     permeability;       // permeability to this field
            int8_t      amplitude;          // intensity of the field pulse generated
            int8_t      changeThreshold;    // threshold to instigate the transition
            int8_t      spawnThreshold;     // threshold to spawn a new stemCell
            uint8_t     fieldType;          // which field the pulse should be of
            uint8_t     nextType;           // type to transition to
            Direction   direction;          // direction in which to spawn a new stemCell
        } gene[fieldsNumber];
    } autosome[stemCellsTypes];
    struct{
        int16_t initialFieldValues[fieldsNumber];         // the initial values of the fields in the first cell
    } allosome;
};

struct StemCell {
    //4*1 bytes = 4 bytes
    uint8_t type;
    int8_t  indices[3];

    //6*8 bytes = 48 bytes
    RelativeStemCellIndex neighbours[6];
    
    //fieldsNumber*2 = 2*fieldsNumber bytes
    int16_t fields[fieldsNumber];
};

struct Embryo{
    uint32_t currentOccupation;
    uint64_t currentSize;
    uint32_t maxStemCells;
    StemCell *stemCells;
    Genome_t genome;
    std::unordered_map<uint32_t, RelativeStemCellIndex> indicesToStemCell;
};

struct Cell{
    uint8_t type;
    int8_t indices[3];
};

struct Body{
    Cell* cells;
    uint32_t size;
};

void generateGenome(Genome_t *genome);

void initializeEmbryo(Embryo *embryo, const Genome_t& genome, uint64_t maxNumber = UINT8_MAX*UINT8_MAX*UINT8_MAX/4);

void reuseEmbryo(Embryo *embryo, const Genome_t& genome);

void copyEmbryo(Embryo* dest, const Embryo* src);

void deleteEmbryo(Embryo *embryo);

RelativeStemCellIndex findStemCellByIndices(Embryo* embryo, int8_t x, int8_t y, int8_t z);

void checkForFieldsSources(Embryo *embryo, RelativeStemCellIndex me);

void diffuse(Embryo *embryo, RelativeStemCellIndex me, RelativeStemCellIndex previous);

void checkForSpeciation(Embryo *embryo, RelativeStemCellIndex me);

void checkForSpawn(Embryo *embryo, RelativeStemCellIndex me);

//Erases the uncoding bits in directions, returns 0 for working genomes and -1 for defective ones.
int sanitizeGenome(Genome_t *genome);

void mutateGenome(Genome_t *genome, float mutationProbability);

//Executes a complete cycle of check for sources-> diffuse -> speciation -> spawn for all StemCells in the embryo
//You shouldn't need to use the individual functions in most cases
void developEmbryo(Embryo* embryo);

void birthBody(const Embryo& embryo, Body* body);

void copyBody(Body* dest, const Body* src);

int geneticDistance(const Genome_t& first, const Genome_t& second);
