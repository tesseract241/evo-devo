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

RelativeCellIndex findCellByIndices(Body* body, int8_t x, int8_t y, int8_t z){
    auto neighbour = body->indicesToCell.find((uint8_t(z)<<16) | (uint8_t(y)<<8) | uint8_t(x));
    if(neighbour!=body->indicesToCell.end()){
        return neighbour->second;
    } else{
        return -1;
    }
}

Cell* newCell(Body *body, uint8_t type, int8_t x, int8_t y, int8_t z){
    body->currentOccupation++;
    if(body->currentOccupation>=body->currentSize){
        Cell *temp = new Cell[2*body->currentSize];
        std::memcpy(temp, body->cells, (body->currentOccupation)*sizeof(Cell));
        delete body->cells;
        body->cells=temp;
        body->currentSize*=2;
    }
    Cell *cell = &(body->cells[body->currentOccupation-1]);
    uint8_t dummy[4] = {type, uint8_t(x), uint8_t(y), uint8_t(z)};
    std::memcpy(&(cell->type), dummy, 4);
    body->indicesToCell[(uint8_t(z)<<16)|(uint8_t(y)<<8)|uint8_t(x)] = body->currentOccupation-1;
    std::memset(cell->neighbours, -1, 6*sizeof(RelativeCellIndex));
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
            RelativeCellIndex neighbour = findCellByIndices(body, neighbourIndices[0], neighbourIndices[1], neighbourIndices[2]);
            if(neighbour!=-1){
                cell->neighbours[i]=neighbour;
                body->cells[neighbour].neighbours[i^1] = cell - body->cells;
            }
        }
    }
    return cell;
}

void generateGenome(Genome_t *genome){
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_int_distribution<uint64_t> dis64(0, UINT64_MAX);
    std::uniform_int_distribution<uint32_t> dis32(0, UINT32_MAX);
    std::uniform_int_distribution<uint8_t> dist(0, cellsTypes-1);
    std::uniform_int_distribution<uint8_t> dis8(1, UINT8_MAX);
    std::uniform_int_distribution<uint8_t> disd(0, 7);
    div_t div = std::div(fieldsNumber, sizeof(uint64_t));
    for(int i=0;i<cellsTypes;++i){
        for(int j=0;j<fieldsNumber;++j){
            uint32_t *ptr = (uint32_t*)(genome->autosome[i].gene + j);
            *ptr = dis32(gen);
            genome->autosome[i].gene[j].nextType = dist(gen);
            genome->autosome[i].gene[j].direction = directions[disd(gen)];
        }
    }
    for(int i=0;i<div.quot;++i){
        uint64_t *ptr = (uint64_t*)(genome->allosome.mass+i*sizeof(uint64_t));
        *ptr = dis64(gen);
    }
    for(int i=0;i<div.rem;++i){
        genome->allosome.mass[div.quot*sizeof(uint64_t)+i] = dis8(gen);
    }
}

void initializeBody(Body *body, const Genome_t& genome){
    body->cells= new Cell[16];
    body->currentOccupation=0;
    body->currentSize=16;
    body->indicesToCell.reserve(16);
    body->genome = genome;
    newCell(body, 0, 0, 0, 0);
}

void reuseBody(Body *body, const Genome_t& genome){
    if(!(body->cells)){
        initializeBody(body, genome);
    } else{
        body->currentOccupation = 0;
        body->genome = genome;
        newCell(body, 0, 0, 0, 0);
    }
}

void copyBody(Body* dest, Body* src){
    if(dest->currentSize < src->currentSize){
        delete dest->cells;
        dest->cells = new Cell[src->currentSize];
    }
    std::memcpy(dest->cells, src->cells, src->currentOccupation * sizeof(Cell));
    dest->currentOccupation = src->currentOccupation;
    dest->currentSize = src->currentSize;
    std::memcpy(&(dest->genome), &(src->genome), sizeof(Genome_t));
    dest->indicesToCell = src->indicesToCell;
}

void deleteBody(Body *body){
    delete body->cells;
    body->currentOccupation=0;
    body->currentSize=0;
    body->indicesToCell.clear();
}

Cell* spawnCell(Body *body, RelativeCellIndex parent, Direction d){
    int16_t indices[3] = {body->cells[parent].indices[0], body->cells[parent].indices[1], body->cells[parent].indices[2]};
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

    Cell *child= newCell(body, body->cells[parent].type, indices[0], indices[1], indices[2]);
    child->neighbours[d^1] = parent;

    for(int i=0;i<fieldsNumber;++i){
        for(int j=0;j<=BACK;++j){
            RelativeCellIndex neighbour = child->neighbours[j];
            if(neighbour!=-1){
                child->fields[i] += body->cells[neighbour].fields[i];
            }
        }
    }
    return child;
}

void pulseField(Body *body, RelativeCellIndex me, int fieldIndex, int16_t intensity, Direction d){
   body->cells[me].fields[fieldIndex]+=intensity;
   if(intensity<2){
        return;
   }
   //The -1 is due to us skipping the cell the pulse is coming from
    uint8_t neighboursCount = -1;
    for(int i=0;i<=BACK;++i){
        if(body->cells[me].neighbours[i]!=-1) ++neighboursCount;
    }
   for(int i=0;i<=BACK;++i){
       if(body->cells[me].neighbours[i]!=-1){
           //This if condition prevents infinite back-propagation
            if(i!=(d^1)){
                pulseField(body, body->cells[me].neighbours[i], fieldIndex, intensity/neighboursCount, (Direction) i);
            }
       }
   }
}

void pulseField(Body *body, RelativeCellIndex me, int fieldIndex, int8_t intensity){
    body->cells[me].fields[fieldIndex] += intensity;
    uint8_t neighboursCount = 0;
    for(int i=0;i<=BACK;++i){
        if(body->cells[me].neighbours[i]!=-1) ++neighboursCount;
    }
    for(int i=0;i<=BACK;++i){
       if(body->cells[me].neighbours[i]!=-1){
            pulseField(body, body->cells[me].neighbours[i], fieldIndex, intensity/neighboursCount, (Direction) i);
        }
    }
}


void checkForSpeciation(Body *body, RelativeCellIndex me){
    for(int i=0;i<fieldsNumber;++i){
        int16_t threshold = (body->genome.autosome[body->cells[me].type].gene[i].changeThreshold) << 8;
        if(threshold==0) continue;
        if(threshold>0){
            if(body->cells[me].fields[i] > threshold){
                body->cells[me].type = body->genome.autosome[body->cells[me].type].gene[i].nextType;
                break;
            }
        }
        else{
            if(body->cells[me].fields[i] < threshold){
                body->cells[me].type = body->genome.autosome[body->cells[me].type].gene[i].nextType;
                break;
            }
        }
    }
}

void checkForPulse(Body *body, RelativeCellIndex me){
    for(int i=0;i<fieldsNumber;++i){
        int8_t amplitude = body->genome.autosome[body->cells[me].type].gene[i].amplitude;
        if(amplitude){
            pulseField(body, me, i, amplitude);
        }
    }
}

void checkForSpawn(Body *body, RelativeCellIndex me){
    for(int i=0;i<fieldsNumber;++i){
        Direction d = body->genome.autosome[body->cells[me].type].gene[i].direction;
        //We only spawn new cells in empty spaces
        if(body->cells[me].neighbours[d]==-1){
            int16_t threshold = (body->genome.autosome[body->cells[me].type].gene[i].spawnThreshold) << 8;
            if(threshold==0) continue;
            if(threshold>0){
                if(body->cells[me].fields[i] > threshold){
                    spawnCell(body, me, d);
                }
            }
            else{
                if(body->cells[me].fields[i] < threshold){
                    spawnCell(body, me, d);
                }
            }
        }
    }
}

void diffuse(Body *body, RelativeCellIndex me){
    for(int i=0;i<fieldsNumber;++i){
        int16_t myField = body->cells[me].fields[i];
        //We save the current state of this cell's permeability so that we can do the exchanges serially
        uint8_t myPermeability = body->genome.autosome[body->cells[me].type].gene[i].permeability;
        if(myPermeability==0){
            return;
        }
        uint8_t mass = body->genome.allosome.mass[i];
        for(uint8_t j=0;j<=BACK;++j){
            RelativeCellIndex neighbour = body->cells[me].neighbours[j];
            if(neighbour!=-1){
                uint8_t theirPermeability = body->genome.autosome[body->cells[neighbour].type].gene[i].permeability;
                if(theirPermeability){
                    //We divide by 4 instead of 2 because diffusion will happen again when called for the neighbour
                    int16_t fieldDelta = (myPermeability*myField - theirPermeability*body->cells[neighbour].fields[i])/(2*myPermeability*theirPermeability*mass);
                    body->cells[me].fields[i]+=fieldDelta;
                    body->cells[neighbour].fields[i]-=fieldDelta;
                }
            }
        }
    }
}

int sanitizeGenome(Genome_t *genome){
    for(int i=0;i<cellsTypes;++i){
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
    for(int i=0;i<cellsTypes;++i){
        for(int j=0;j<fieldsNumber;++j){
            //Checks if nextType is one of the ones available
            if(genome->autosome[i].gene[j].nextType>cellsTypes-1){
                return -1;
            }
        }
    }
    return 0;
}

void mutateGenome(Genome_t *genome, float mutationProbability){
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<float> prob(0, 1);
    std::uniform_int_distribution<uint8_t> dist(0, cellsTypes-1);
    std::uniform_int_distribution<uint8_t> dis8(0, 7);
    std::uniform_int_distribution<uint8_t> dism(1, 7);
    std::uniform_int_distribution<uint8_t> disb(0, 1);
    for(int i=0;i<cellsTypes;++i){
        for(int j=0;j<fieldsNumber;++j){
            uint8_t *ptr = (uint8_t*) &(genome->autosome[i].gene[j]);
            for(int k=0;k<4;++k){
                if(prob(gen)<mutationProbability){
                    ptr[k]^=(1<<dis8(gen));
                }
            }
            if(prob(gen)<mutationProbability){
                genome->autosome[i].gene[j].nextType = dist(gen);
            }
            if(prob(gen)<mutationProbability){
                uint8_t currentDirection = genome->autosome[i].gene[j].direction;
                if(currentDirection>7){
                    currentDirection-=2;
                }
                genome->autosome[i].gene[j].direction = directions[(currentDirection+2*disb(gen)-1)%8];
            }
        }
    }
    for(int i=0;i<fieldsNumber;++i){
        if(prob(gen)<mutationProbability){
            genome->allosome.mass[i]^=(1<<dism(gen));
        }
    }
}

void developBody(Body* body){
    for(int i=0;i<body->currentOccupation;++i){
        checkForPulse(body, i);
    }
    for(int i=0;i<body->currentOccupation;++i){
        diffuse(body, i);
    }
    for(int i=0;i<body->currentOccupation;++i){
        checkForSpeciation(body, i);
    }
    for(int i=0;i<body->currentOccupation;++i){
        checkForSpawn(body, i);
    }
}

int geneticDistance(const Genome_t& first, const Genome_t& second){
    int sum = 0;
    for(int i=0;i<cellsTypes;++i){
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
        sum+=(first.allosome.mass[j] - second.allosome.mass[j])*(first.allosome.mass[j] - second.allosome.mass[j]);
    }
    return sum;
}
