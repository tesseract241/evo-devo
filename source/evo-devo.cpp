#include <cassert>
#include <cmath>
#include <cstring>
#include <evo-devo.hpp>

Cell* findCellByIndices(Body* body, int8_t x, int8_t y, int8_t z){
    auto neighbour = body->indicesToCell.find((uint8_t(z)<<16) | (uint8_t(y)<<8) | uint8_t(x));
    if(neighbour!=body->indicesToCell.end()){
        return body->cells + neighbour->second;
    } else{
        return NULL;
    }
}

Cell* newCell(Body *body, uint8_t type, int8_t x, int8_t y, int8_t z){
    body->currentOccupation++;
    if(body->currentOccupation>=body->currentSize){
        Cell *temp = new Cell[2*body->currentSize];
        std::memcpy(temp, body->cells, (body->currentOccupation)*sizeof(Cell));
        delete body->cells;
        body->cells=temp;
    }
    Cell *cell = &(body->cells[body->currentOccupation]);
    uint8_t dummy[4] = {type, uint8_t(x), uint8_t(y), uint8_t(z)};
    std::memcpy(&(cell->type), dummy, 4);
    body->indicesToCell[(uint8_t(z)<<16)|(uint8_t(y)<<8)|uint8_t(x)] = body->currentOccupation;
    return cell;
}

void initializeBody(Body *body, const Genome_t& genome){
    body->cells= new Cell[16];
    body->currentOccupation=-1;
    body->currentSize=16;
    body->indicesToCell.reserve(16);
    body->genome = genome;
    newCell(body, 0, 0, 0, 0);
}

void reuseBody(Body *body, const Genome_t& genome){
    if(!(body->cells)){
        initializeBody(body, genome);
    } else{
        body->currentOccupation = -1;
        body->genome = genome;
        newCell(body, 0, 0, 0, 0);
    }
}

void deleteBody(Body *body){
    delete body->cells;
    body->currentOccupation=0;
    body->currentSize=0;
    body->indicesToCell.clear();
}

Cell* spawnCell(Body *body, Cell *parent, Direction d){
    int16_t indices[3] = {parent->indices[0],parent->indices[1],parent->indices[2]};
    //If the direction is relative (INWARDS/OUTWARDS) we need to describe it as an absolute direction
    //by looking at which coordinate is the greatest in absolute value
    if(d&6){
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

    Cell *child= newCell(body, parent->type, indices[0], indices[1], indices[2]);
    child->type = parent->type;
    *(&(child->up)+(d^1)) = parent;
    for(int i=0;i<=BACK;++i){
        int16_t neighbourIndices[3] = {indices[0], indices[1], indices[2]};
        if(i!=(d^1)){
            neighbourIndices[(i&6)>>1]-= (i&1)*2 - 1;
            bool skip = false;
            for(int j=0;j<3;++j){
                if(neighbourIndices[i]<INT8_MIN || neighbourIndices[i]> INT8_MAX){
                    skip = true;
                    break;
                }
            }
            if(!skip){
                *(&(child->up)+i)=findCellByIndices(body, neighbourIndices[0], neighbourIndices[1], neighbourIndices[2]);
            }
        }
    }

    for(int i=0;i<fieldsNumber;++i){
        for(int j=0;j<=BACK;++j){
            Cell* neighbour = *(&(child->up)+j);
            if(neighbour){
                child->fields[i] += neighbour->fields[i];
            }
        }
    }
    return child;
}

void pulseField(Cell* me, int fieldIndex, int16_t intensity, Direction d){
   me->fields[fieldIndex]+=intensity;
   if(intensity==1){
        return;
   }
   //The -1 is due to us skipping the cell the pulse is coming from
    uint8_t neighboursCount = -1;
    for(int i=0;i<=BACK;++i){
        if(*(&(me->up)+i)) ++neighboursCount;
    }
   for(int i=0;i<=BACK;++i){
       if(*(&(me->up)+i)){
           //This if condition prevents infinite back-propagation
            if(i!=(d^1)){
                pulseField(me->up+i, fieldIndex, intensity/neighboursCount, (Direction) i);
            }
       }
   }
}

void pulseField(Cell* me, int fieldIndex, int8_t intensity){
    me->fields[fieldIndex] += intensity;
    uint8_t neighboursCount = 0;
    for(int i=0;i<=BACK;++i){
        if(*(&(me->up)+i)) ++neighboursCount;
    }
    for(int i=0;i<=BACK;++i){
       if(*(&(me->up)+i)){
            pulseField(me->up+i, fieldIndex, intensity/neighboursCount, (Direction) i);
        }
    }
}


void checkForSpeciation(Body *body, Cell* me){
    for(int i=0;i<fieldsNumber;++i){
        int16_t threshold = (body->genome.autosome[me->type].gene[i].changeThreshold) << 8;
        if(threshold==0) continue;
        if(threshold>0){
            if(me->fields[i] > threshold){
                me->type = body->genome.autosome[me->type].gene[i].nextType;
                break;
            }
        }
        else{
            if(me->fields[i] < threshold){
                me->type = body->genome.autosome[me->type].gene[i].nextType;
                break;
            }
        }
    }
}

void checkForPulse(Body *body, Cell* me){
    for(int i=0;i<fieldsNumber;++i){
        int8_t amplitude = body->genome.autosome[me->type].gene[i].amplitude;
        if(amplitude){
            pulseField(me, i, amplitude);
        }
    }
}

void checkForSpawn(Body *body, Cell* me){
    for(int i=0;i<fieldsNumber;++i){
        Direction d = body->genome.autosome[me->type].gene[i].direction;
        //We only spawn new cells in empty spaces
        if(!*(&(me->up)+d)){
            int16_t threshold = (body->genome.autosome[me->type].gene[i].spawnThreshold) << 8;
            if(threshold==0) continue;
            if(threshold>0){
                if(me->fields[i] > threshold){
                    spawnCell(body, me, d);
                }
            }
            else{
                if(me->fields[i] < threshold){
                    spawnCell(body, me, d);
                }
            }
        }
    }
}

void diffuse(Body *body, Cell *me){
    for(int i=0;i<fieldsNumber;++i){
        int16_t myField = me->fields[i];
        //We save the current state of this cell's permeability so that we can do the exchanges serially
        uint8_t myPermeability = body->genome.autosome[me->type].gene[i].permeability;
        uint8_t mass = body->genome.allosome.mass[i];
        for(uint8_t j=0;j<=BACK;++j){
            Cell* neighbour = *(&(me->up)+j);
            if(neighbour){
                uint8_t theirPermeability = body->genome.autosome[neighbour->type].gene[i].permeability;
                //We divide by 4 instead of 2 because diffusion will happen again when called for the neighbour
                int16_t fieldDelta = (myPermeability*myField - theirPermeability*neighbour->fields[i])/(2*myPermeability*theirPermeability*mass);
                    me->fields[i]+=fieldDelta;
                    neighbour->fields[i]-=fieldDelta;
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

void developBody(Body* body){
    for(int i=0;i<body->currentOccupation;++i){
        checkForPulse(body, body->cells+i);
    }
    for(int i=0;i<body->currentOccupation;++i){
        diffuse(body, body->cells+i);
    }
    for(int i=0;i<body->currentOccupation;++i){
        checkForSpeciation(body, body->cells+i);
    }
    for(int i=0;i<body->currentOccupation;++i){
        checkForSpawn(body, body->cells+i);
    }
}
