#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <random>
#include <evo-devo.hpp>

const Direction directions[8] = {
         RIGHT,
         LEFT,
         UP,
         DOWN,
         FORWARD,
         BACK,
         OUTWARDS,
         INWARDS
};

RelativeStemCellIndex findStemCellByIndices(Embryo* embryo, int8_t x, int8_t y, int8_t z){
    auto neighbour = embryo->indicesToStemCell.find((uint8_t(z)<<16) | (uint8_t(y)<<8) | uint8_t(x));
    if(neighbour!=embryo->indicesToStemCell.end()){
        return neighbour->second;
    } else{
        return -1;
    }
}

StemCell* newStemCell(Embryo *embryo, uint8_t type, int8_t x, int8_t y, int8_t z){
    embryo->currentOccupation++;
    if(embryo->currentOccupation>=embryo->currentSize){
        StemCell *temp = new StemCell[2*embryo->currentSize];
        std::memcpy(temp, embryo->stemCells, (embryo->currentOccupation)*sizeof(StemCell));
        delete embryo->stemCells;
        embryo->stemCells=temp;
        embryo->currentSize*=2;
    }
    StemCell *stemCell = &(embryo->stemCells[embryo->currentOccupation-1]);
    uint8_t dummy[4] = {type, uint8_t(x), uint8_t(y), uint8_t(z)};
    std::memcpy(&(stemCell->type), dummy, 4);
    embryo->indicesToStemCell[(uint8_t(z)<<16)|(uint8_t(y)<<8)|uint8_t(x)] = embryo->currentOccupation-1;
    std::memset(stemCell->neighbours, -1, 6*sizeof(RelativeStemCellIndex));
    for(int i=0;i<=BACK;++i){
        int16_t neighbourIndices[3] = {x, y, z};
        neighbourIndices[(i&6)>>1]-= (i&1)*2 - 1;
        bool skip = false;
        for(int j=0;j<3;++j){
            if(neighbourIndices[j]<INT8_MIN || neighbourIndices[j]> INT8_MAX){
                skip = true;
                break;
            }
        }
        if(!skip){
            RelativeStemCellIndex neighbour = findStemCellByIndices(embryo, neighbourIndices[0], neighbourIndices[1], neighbourIndices[2]);
            if(neighbour!=-1){
                stemCell->neighbours[i]=neighbour;
                embryo->stemCells[neighbour].neighbours[i^1] = stemCell - embryo->stemCells;
            }
        }
    }
    return stemCell;
}

void generateGenome(Genome_t *genome){
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_int_distribution<uint64_t> dis64(0,        UINT64_MAX);
    std::uniform_int_distribution<uint32_t> dis32(0,        UINT32_MAX);
    std::uniform_int_distribution<uint8_t>  dist( 0,        stemCellsTypes-1);
    std::uniform_int_distribution<uint8_t>  disf( 0,        fieldsNumber-1);
    std::uniform_int_distribution<uint8_t>  dis8( 1,        UINT8_MAX);
    std::uniform_int_distribution<int8_t>   disi( INT8_MIN, INT8_MAX);
    std::uniform_int_distribution<uint8_t>  disd( 0,        7);
    div_t div = std::div((int)(fieldsNumber * sizeof(int16_t)), sizeof(uint64_t));
    for(int i=0;i<stemCellsTypes;++i){
        for(int j=0;j<fieldsNumber;++j){
            uint32_t *ptr = (uint32_t*)(genome->autosome[i].gene + j);
            *ptr = dis32(gen);
            genome->autosome[i].gene[j].spawnThreshold  = disi(gen);
            genome->autosome[i].gene[j].fieldType       = disf(gen);
            genome->autosome[i].gene[j].nextType        = dist(gen);
            genome->autosome[i].gene[j].direction       = directions[disd(gen)];
        }
    }
    for(int i=0;i<div.quot;++i){
        uint64_t *ptr = (uint64_t*)(genome->allosome.initialFieldValues+i*sizeof(uint64_t));
        *ptr = dis64(gen);
    }
    for(int i=0;i<div.rem;++i){
        genome->allosome.initialFieldValues[div.quot*sizeof(uint64_t) + i] = dis8(gen);
    }
}

void initializeEmbryo(Embryo *embryo, const Genome_t& genome, uint64_t maxNumber){
    if(maxNumber!=0){
        embryo->maxStemCells = maxNumber;
        embryo->stemCells= new StemCell[16];
        embryo->currentOccupation=0;
        embryo->currentSize=16;
        embryo->indicesToStemCell.reserve(16);
        embryo->genome = genome;
        newStemCell(embryo, 0, 0, 0, 0);
        std::memcpy(embryo->stemCells[0].fields, genome.allosome.initialFieldValues, fieldsNumber*sizeof(int16_t));
    }
}

void reuseEmbryo(Embryo *embryo, const Genome_t& genome){
    if(!(embryo->stemCells)){
        initializeEmbryo(embryo, genome);
    } else{
        embryo->currentOccupation = 0;
        embryo->genome = genome;
        newStemCell(embryo, 0, 0, 0, 0);
        std::memcpy(embryo->stemCells[0].fields, genome.allosome.initialFieldValues, fieldsNumber*sizeof(int16_t));
    }
}

void copyEmbryo(Embryo* dest, const Embryo* src){
    if(dest->currentSize < src->currentSize){
        delete dest->stemCells;
        dest->stemCells = new StemCell[src->currentSize];
    }
    std::memcpy(dest->stemCells, src->stemCells, src->currentOccupation * sizeof(StemCell));
    dest->currentOccupation = src->currentOccupation;
    dest->currentSize = src->currentSize;
    std::memcpy(&(dest->genome), &(src->genome), sizeof(Genome_t));
    dest->indicesToStemCell = src->indicesToStemCell;
}

void deleteEmbryo(Embryo *embryo){
    delete embryo->stemCells;
    embryo->currentOccupation=0;
    embryo->currentSize=0;
    embryo->indicesToStemCell.clear();
}

StemCell* spawnStemCell(Embryo *embryo, RelativeStemCellIndex parent, Direction d){
    int16_t indices[3] = {embryo->stemCells[parent].indices[0], embryo->stemCells[parent].indices[1], embryo->stemCells[parent].indices[2]};
    //If the direction is relative (INWARDS/OUTWARDS) we need to describe it as an absolute direction
    //by looking at which coordinate is the greatest in absolute value
    if(d&OUTWARDS){
        uint8_t max = 0;
        Direction newDirection;
        for(int i=0;i<3;++i){
            if(std::abs(indices[i]>max)){
                max = std::abs(indices[i]);
                newDirection = (Direction) ((1<<i)  &  (6 | (uint8_t(indices[i])>>7)));
//                                                         sign of index moved to the LSB
            }
        }
        d = newDirection;
    }
    //This adds/subtracts 1 from the correct index based on direction to get the child's coordinates
    indices[(d&6)>>1]-= (d&1)*2 - 1;
    for(int i=0;i<3;++i){
        if(indices[i]<INT8_MIN || indices[i]>INT8_MAX){
            return NULL;
        }
    }

    StemCell *child= newStemCell(embryo, embryo->stemCells[parent].type, indices[0], indices[1], indices[2]);
    child->neighbours[d^1] = parent;

    for(int i=0;i<fieldsNumber;++i){
        for(int j=0;j<=BACK;++j){
            RelativeStemCellIndex neighbour = child->neighbours[j];
            if(neighbour!=-1){
                child->fields[i] += embryo->stemCells[neighbour].fields[i];
            }
        }
    }
    return child;
}

void checkForSpeciation(Embryo *embryo, RelativeStemCellIndex me){
    for(int i=0;i<fieldsNumber;++i){
        int16_t threshold = (embryo->genome.autosome[embryo->stemCells[me].type].gene[i].changeThreshold) << 8;
        if(threshold==0) continue;
        if(threshold>0){
            if(embryo->stemCells[me].fields[i] > threshold){
                embryo->stemCells[me].type = embryo->genome.autosome[embryo->stemCells[me].type].gene[i].nextType;
                break;
            }
        }
        else{
            if(embryo->stemCells[me].fields[i] < threshold){
                embryo->stemCells[me].type = embryo->genome.autosome[embryo->stemCells[me].type].gene[i].nextType;
                break;
            }
        }
    }
}

void checkForFieldsSources(Embryo *embryo, RelativeStemCellIndex me){
    for(int i=0;i<fieldsNumber;++i){
        int16_t pulseThreshold = embryo->genome.autosome[embryo->stemCells[me].type].gene[i].pulseThreshold<<8;
        uint8_t fieldType = embryo->genome.autosome[embryo->stemCells[me].type].gene[i].fieldType;
        if(pulseThreshold==0) continue;
        if(
                (pulseThreshold > 0 && embryo->stemCells[me].fields[i] > pulseThreshold) ||
                (pulseThreshold < 0 && embryo->stemCells[me].fields[i] < pulseThreshold) 
        //Note that this will create a new element if it doesn't exist
        ){
            embryo->stemCells[me].fields[fieldType]+=embryo->genome.autosome[embryo->stemCells[me].type].gene[i].amplitude;
            //embryo->sources[me] |= 1<<fieldType;
        }
    }
}

void checkForSpawn(Embryo *embryo, RelativeStemCellIndex me){
    for(int i=0;i<fieldsNumber;++i){
        Direction d = embryo->genome.autosome[embryo->stemCells[me].type].gene[i].direction;
        //We only spawn new stemCells in empty spaces
        if(embryo->stemCells[me].neighbours[d]==-1){
            int16_t threshold = (embryo->genome.autosome[embryo->stemCells[me].type].gene[i].spawnThreshold) << 8;
            if(threshold==0) continue;
            if(threshold>0){
                if(embryo->stemCells[me].fields[i] > threshold){
                    spawnStemCell(embryo, me, d);
                }
            }
            else{
                if(embryo->stemCells[me].fields[i] < threshold){
                    spawnStemCell(embryo, me, d);
                }
            }
        }
    }
}

//void diffuse(Embryo *embryo, RelativeStemCellIndex me, int field, uint8_t range, RelativeStemCellIndex previous){
//    uint8_t myPermeability = embryo->genome.autosome[embryo->stemCells[me].type].gene[field].permeability;
//    if(myPermeability==0){
//        return;
//    }
//    int16_t myField = embryo->stemCells[me].fields[field];
//    //We save the current state of this stemCell's field so that we can do the exchanges serially
//    for(uint8_t j=0;j<=BACK;++j){
//        RelativeStemCellIndex neighbour = embryo->stemCells[me].neighbours[j];
//        if(neighbour==-1) continue;
//        uint8_t theirPermeability = embryo->genome.autosome[embryo->stemCells[neighbour].type].gene[field].permeability;
//        //We check that the permeability is not zero, that the neighbour is not the stemCell we come from and not in the list of sources for this field
//        if(theirPermeability && neighbour!=previous){
//            auto neighbourIterator = embryo->sources.find(neighbour);
//            if(neighbourIterator!=embryo->sources.end() && (neighbourIterator->second&(1<<field))){
//                int16_t fieldDelta = (myPermeability*myField - theirPermeability*embryo->stemCells[neighbour].fields[field])/(2*(myPermeability+theirPermeability));
//                embryo->stemCells[me].fields[field]+=fieldDelta;
//                embryo->stemCells[neighbour].fields[field]-=fieldDelta;
//                if(range!=1) diffuse(embryo, neighbour, field, range-1, me);
//            }
//        }
//    }
//}

void diffuse(Embryo *embryo, RelativeStemCellIndex me, RelativeStemCellIndex previous){
    for(int i=0;i<fieldsNumber;++i){
        uint8_t myPermeability = embryo->genome.autosome[embryo->stemCells[me].type].gene[i].permeability;
        if(myPermeability==0){
            return;
        }
        int16_t myField = embryo->stemCells[me].fields[i];
        for(uint8_t j=0;j<=BACK;++j){
            RelativeStemCellIndex neighbour = embryo->stemCells[me].neighbours[j];
            if(neighbour==-1) continue;
            uint8_t theirPermeability = embryo->genome.autosome[embryo->stemCells[neighbour].type].gene[i].permeability;
            //We check that the permeability is not zero and that the neighbour is not the stemCell we come from
            if(theirPermeability && neighbour!=previous){
                //We divide by two because neighbour will repeat the exchange
                int16_t fieldDelta = (myPermeability*myField - theirPermeability*embryo->stemCells[neighbour].fields[i])/(2*(myPermeability+theirPermeability));
                embryo->stemCells[me].fields[i]+=fieldDelta;
                embryo->stemCells[neighbour].fields[i]-=fieldDelta;
            }
        }
    }
}

int sanitizeGenome(Genome_t *genome){
    for(int i=0;i<stemCellsTypes;++i){
        for(int j=0;j<fieldsNumber;++j){
            //Turns off the non-coding bits
            genome->autosome[i].gene[j].direction = (Direction) (genome->autosome[i].gene[j].direction&15);
            Direction d = genome->autosome[i].gene[j].direction;
            //Quick and dirty way to count the number of ON bits, excluding the LSB
            //The if condition fails if d isn't a valid Direction
            if((((d&2)>>1) + ((d&4)>>2) + ((d&8)>>3))>1){
                return -1;
            }
        }
    }
    for(int i=0;i<stemCellsTypes;++i){
        for(int j=0;j<fieldsNumber;++j){
            //Checks if nextType is one of the ones available
            if(genome->autosome[i].gene[j].nextType>stemCellsTypes-1){
                return -1;
            }
        }
    }
    return 0;
}

void mutateGenome(Genome_t *genome, float mutationProbability){
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<float>   prob (0, 1);
    std::uniform_int_distribution<uint8_t>  dist (0, stemCellsTypes-1);
    std::uniform_int_distribution<uint8_t>  disft(0, fieldsNumber-1);
    std::uniform_int_distribution<uint8_t>  dis8 (0, 7);
    std::uniform_int_distribution<uint8_t>  disf (0, 15);
    std::uniform_int_distribution<uint8_t>  disb (0, 1);
    for(int i=0;i<stemCellsTypes;++i){
        for(int j=0;j<fieldsNumber;++j){
            uint8_t *ptr = (uint8_t*) &(genome->autosome[i].gene[j]);
            for(int k=0;k<5;++k){
                if(prob(gen)<mutationProbability){
                    ptr[k]^=(1<<dis8(gen));
                }
            }
            if(prob(gen)<mutationProbability){
                genome->autosome[i].gene[j].fieldType= disft(gen);
            }
            if(prob(gen)<mutationProbability){
                genome->autosome[i].gene[j].nextType = dist(gen);
            }
            if(prob(gen)<mutationProbability){
                uint8_t currentDirection = genome->autosome[i].gene[j].direction;
                //We want the direction to move one step from where it is, either forward or backward depending on the binary random number
                //To do this, we use the directions array and the modulus operation. Notice that when going out of bounds this cycles back to the other side
                if(currentDirection>7){
                    currentDirection-=2;
                }
                genome->autosome[i].gene[j].direction = directions[(currentDirection+2*disb(gen)+7)%8];
            }
        }
    }
    for(int i=0;i<fieldsNumber;++i){
        if(prob(gen)<mutationProbability){
            genome->allosome.initialFieldValues[i]^=(1<<disf(gen));
        }
    }
}

void developEmbryo(Embryo* embryo, uint8_t range){
    for(uint64_t i=0;i<embryo->currentOccupation;++i){
        checkForFieldsSources(embryo, i);
    }
        for(uint64_t i=0;i<embryo->currentOccupation;++i){
            diffuse(embryo, i, -1);
        }
    for(uint64_t i=0;i<embryo->currentOccupation;++i){
        checkForSpeciation(embryo, i);
    }
    for(uint64_t i=0;i<embryo->currentOccupation;++i){
        if(embryo->currentOccupation<embryo->maxStemCells){
            checkForSpawn(embryo, i);
        }
    }
}

void birthBody(const Embryo& embryo, Body* body){
    if(!body->cells){
        body->cells = new Cell[embryo.currentOccupation];
        body->size = embryo.currentOccupation;
    } else {
        if(body->size < embryo.currentOccupation){
            delete body->cells;
            body->cells = new Cell[embryo.currentOccupation];
            body->size = embryo.currentOccupation;
        }
    }
    for(uint64_t i=0;i<embryo.currentOccupation;++i){
        std::memcpy(body->cells + i, embryo.stemCells + i, 4);
    }
    body->occupation = embryo.currentOccupation;
}

void copyBody(Body* dest, const Body* src){
    if(dest->size < src->occupation){
        delete dest->cells;
        dest->cells = new Cell[src->occupation];
    }
    std::memcpy(dest->cells, src->cells, src->occupation * sizeof(Cell));
    dest->occupation = src->occupation;
    dest->size = dest->occupation;
}

int geneticDistance(const Genome_t& first, const Genome_t& second){
    int sum = 0;
    for(int i=0;i<stemCellsTypes;++i){
        for(int j=0;j<fieldsNumber;++j){
            sum+=(first.autosome[i].gene[j].permeability - second.autosome[i].gene[j].permeability)*(first.autosome[i].gene[j].permeability - second.autosome[i].gene[j].permeability);
            sum+=(first.autosome[i].gene[j].amplitude - second.autosome[i].gene[j].amplitude)*(first.autosome[i].gene[j].amplitude - second.autosome[i].gene[j].amplitude);
            sum+=(first.autosome[i].gene[j].changeThreshold - second.autosome[i].gene[j].changeThreshold)*(first.autosome[i].gene[j].changeThreshold - second.autosome[i].gene[j].changeThreshold);
            sum+=(first.autosome[i].gene[j].spawnThreshold - second.autosome[i].gene[j].spawnThreshold)*(first.autosome[i].gene[j].spawnThreshold - second.autosome[i].gene[j].spawnThreshold);
            sum+=(first.autosome[i].gene[j].nextType != second.autosome[i].gene[j].nextType) ? UINT8_MAX*UINT8_MAX : 0;
            sum+=(first.autosome[i].gene[j].direction != second.autosome[i].gene[j].direction) ? UINT8_MAX*UINT8_MAX : 0;
        }
    }
    for(int j=0;j<fieldsNumber;++j){
        sum+=(first.allosome.initialFieldValues[j] - second.allosome.initialFieldValues[j])*(first.allosome.initialFieldValues[j] - second.allosome.initialFieldValues[j]);
    }
    return sum;
}
