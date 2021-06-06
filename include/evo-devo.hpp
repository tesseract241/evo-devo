#pragma once
#include <cstdint>
#include<unordered_map>

const int fieldsNumber  = 14;   //We want a multiple of 4 + 2 to make the Cell struct tightly packed
const int cellsTypes    = 32;

typedef int64_t RelativeCellIndex; 

enum Direction{
    RIGHT       = 0,    //0000
    LEFT        = 1,    //0001
    UP          = 2,    //0010
    DOWN        = 3,    //0011
    FORWARD     = 4,    //0100
    BACK        = 5,    //0101
    OUTWARDS    = 8,    //1000
    INWARDS     = 9     //1001
};


//The genome completely defines the behaviour of all cells during development
struct Genome_t{
    //Each autosome contains the behaviour of a different cell type
    struct {
        //Each gene contains the amplitude that the cell generates for its field, which type it will transition to if exposed to this field, and the threshold above/below which the transition happens, the direction to spawn a new cell, and the threshold associated with it
        struct {
            uint8_t     permeability;       // permeability to this field
            int8_t      amplitude;          // intensity of the field pulse generated
            int8_t      changeThreshold;    // threshold to instigate the transition
            int8_t      spawnThreshold;     // threshold to spawn a new cell
            uint8_t     nextType;           // type to transition to
            Direction   direction;          // direction in which to spawn a new cell
        } gene[fieldsNumber];
    } autosome[cellsTypes];
    struct{
        uint8_t mass[fieldsNumber];         // mass of the field, determines its dynamics
    } allosome;
};

struct Cell {
    //6*8 bytes = 48 bytes
    RelativeCellIndex neighbours[6];
    
    //4*1 bytes = 4 bytes
    uint8_t type;
    int8_t  indices[3];
    
    //fieldsNumber*2 = 2*fieldsNumber bytes
    int16_t fields[fieldsNumber];
};

struct Body{
    int currentOccupation;
    int currentSize;
    Cell *cells;
    Genome_t genome;
    std::unordered_map<uint32_t, uint32_t> indicesToCell;
};

void generateGenome(Genome_t *genome);

void initializeBody(Body *body, const Genome_t& genome);

void reuseBody(Body *body, const Genome_t& genome);

void copyBody(Body* dest, Body* src);

void deleteBody(Body *body);

RelativeCellIndex findCellByIndices(Body* body, int8_t x, int8_t y, int8_t z);

void checkForFieldsSources(Body *body, RelativeCellIndex me);

void diffuse(Body *body, RelativeCellIndex me);

void checkForSpeciation(Body *body, RelativeCellIndex me);

void checkForSpawn(Body *body, RelativeCellIndex me);

//Erases the uncoding bits in directions, returns 0 for working genomes and -1 for defective ones.
int sanitizeGenome(Genome_t *genome);

void mutateGenome(Genome_t *genome, float mutationProbability);

//Executes a complete cycle of pulse field -> diffuse -> speciation -> spawn for all cells in the body
//You shouldn't need to use the individual functions in most cases
void developBody(Body* body);

int geneticDistance(const Genome_t& first, const Genome_t& second);
