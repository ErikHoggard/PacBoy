--
-- Created by IntelliJ IDEA.
-- User: Erik
-- Date: 4/5/2018
-- Time: 8:13 PM
-- To change this template use File | Settings | File Templates.
--

Filename = "start.state"
ButtonNames = {
    "Up",
    "Down",
    "Left",
    "Right",
}

inputRadius = 6; --input is a grid of tiles centered on pacman.
inputCount = ((2 * inputRadius + 1) ^ 2) + 1
outputCount = #ButtonNames --outputs are the buttons of the NES controller

--pool = {}

Population = 200
DeltaDisjoint = 2.0
DeltaWeights = 0.4
DeltaThreshold = 1.0

StaleSpecies = 15

MutateConnectionsChance = 0.25
PerturbChance = 0.90
CrossoverChance = 0.75
LinkMutationChance = 2.0
NodeMutationChance = 0.50
BiasMutationChance = 0.40
StepSize = 0.1
DisableMutationChance = 0.4
EnableMutationChance = 0.2

TimeoutConstant = 20

MaxNodes = 1000000

print(memory.getmemorydomainlist())
memory.usememorydomain("CIRAM (nametables)")

ghostNear = false

score = 0;
bestScore = 0;
lifeCount = 0;
pelletCount = 0;

pacmanX = 0;
pacmanY = 0;
pacmanTile = 0;

blinkyX = 0;
blinkyY = 0;
blinkyTile = 0;
blinkySprite = 0;

pinkyX = 0;
pinkyY = 0;
pinkyTile = 0;
pinkySprite = 0;

inkyX = 0;
inkyY = 0;
inkyTile = 0;
inkySprite = 0;

clydeX = 0;
clydeY = 0;
clydeTile = 0;
clydeSprite = 0;


function updateCharacterPositions()
    --    only check if Y positions are divisible by 8 to preserve vertical tile position in memory

    --    if (pacmanY == 16 or (pacmanY - 16) % 8 == 0) then
    memory.usememorydomain("RAM")

    pelletCount = memory.readbyte(0x006A)
    lifeCount = memory.readbyte(0x0067)

    pacmanX = memory.readbyte(0x001A)
    pacmanY = memory.readbyte(0x001C)
    if ((pacmanY - 16) % 8 == 0) then
        memory.usememorydomain("RAM")
        pacmanTile = getTileAtCoords(pacmanX, pacmanY)
    end

    memory.usememorydomain("RAM")
    blinkyX = memory.readbyte(0x001E)
    blinkyY = memory.readbyte(0x0020)
    blinkySprite = memory.readbyte(0x0033)
    if ((blinkyY - 16) % 8 == 0) then
        blinkyTile = getTileAtCoords(blinkyX, blinkyY)
    end

    memory.usememorydomain("RAM")
    pinkyX = memory.readbyte(0x0022)
    pinkyY = memory.readbyte(0x0024)
    pinkySprite = memory.readbyte(0x0034)
    if ((pinkyY - 16) % 8 == 0) then
        pinkyTile = getTileAtCoords(pinkyX, pinkyY)
    end

    memory.usememorydomain("RAM")
    inkyX = memory.readbyte(0x0026)
    inkyY = memory.readbyte(0x0028)
    inkySprite = memory.readbyte(0x0035)
    if ((inkyY - 16) % 8 == 0) then
        inkyTile = getTileAtCoords(inkyX, inkyY)
    end

    memory.usememorydomain("RAM")
    clydeX = memory.readbyte(0x002A)
    clydeY = memory.readbyte(0x002C)
    clydeSprite = memory.readbyte(0x0036)
    if ((clydeY - 16) % 8 == 0) then
        clydeTile = getTileAtCoords(clydeX, clydeY)
    end
end

-- returns memory address of a tile at coordinates
function getTileAtCoords(x, y)

    if (y == 16 or (y - 16) % 8 == 0) then

        --        memory.usememorydomain("CIRAM (nametables)")

        local base = 98 -- The top-leftmost walkable tile has a memory address of 98
        local tileSize = 8 -- 8 distance units between tiles, horizontally and vertically
        local rowOffsetMult = 32 -- There are 32 bytes between rows of tiles in memory
        local xOffset = 24 -- The top-leftmost walkable position of pacman/ghosts has an x position of 24
        local yOffset = 16 -- The top-leftmost walkable position of pacman/ghosts has a  y position of 16

        local row = base + ((y - yOffset) / tileSize) * rowOffsetMult
        local column = (x - xOffset) / tileSize
        --        print(column)
        --        print(row)

        local memAddress = math.floor(row + column)
        --        print(memAddress)
        local tile = memory.readbyte(memAddress)
        --        previousTile = tile
        --        print(row)

        return memAddress
        --        return tile
        --    else
        --        return previousTile
    end
end

function getTileValue(address)
    memory.usememorydomain("CIRAM (nametables)")
    tile = memory.readbyte(address)
    return tile
end

function genInputs()
    updateCharacterPositions()

    local inputList = {}

    --    for i = -20, 20, 1 do
    --        for j = -16, 16, 1 do
    for i = -6, 6, 1 do
        for j = -6, 6, 1 do
            inputList[#inputList + 1] = 0
            local x = pacmanX + (j * 8)
            local y = pacmanY + (i * 8)

            if (x >= 24 and x <= 168 and y >= 16 and y <= 208) then
                if ((x - 24) % 8 == 0 and (y - 16) % 8 == 0) then
                    --                memory.usememorydomain("CIRAM (nametables)")
                    local addr = getTileAtCoords(x, y)
                    local tile = getTileValue(addr)
                    --                local tile = memory.readbyte(addr)
                    -- detect empty walkable space
                    if (tile == 7 or tile == 0) then
                        inputList[#inputList] = 0.75
                    end
                    -- detect dots
                    if (tile == 3) then
                        inputList[#inputList] = 1
                    end
                    if (tile == 1 or tile == 2) then
                        inputList[#inputList] = 1
                    end
                    -- detect ghosts
                    ghostNear = false
                    if ((addr == blinkyTile and (blinkySprite ~= 30 and blinkySprite ~= 31)) or
                            (addr == pinkyTile and (pinkySprite ~= 30 and pinkySprite ~= 31)) or
                            (addr == inkyTile and (inkySprite ~= 30 and inkySprite ~= 31)) or
                            (addr == clydeTile and (clydeSprite ~= 30 and clydeSprite ~= 31))) then
                        inputList[#inputList] = -1

                    else
                        if ((addr == blinkyTile and (blinkySprite == 30 and blinkySprite == 31)) or
                                (addr == pinkyTile and (pinkySprite ~= 30 and pinkySprite == 31)) or
                                (addr == inkyTile and (inkySprite == 30 and inkySprite == 31)) or
                                (addr == clydeTile and (clydeSprite == 30 and clydeSprite == 31))) then
                            inputList[#inputList] = 1
                        end
                    end
                end
            end
        end
    end

    return inputList
end


function sigmoid(x)
    return 2 / (1 + math.exp(-4.9 * x)) - 1
end

function sigFunction(x)
    sigmoid = 1 / (1 + math.exp(-x))
    --    return 2 / (1 + math.exp(-4.9 * x)) - 1
end

-- Increment the current innovation number
function newInnovation()
    pool.innovation = pool.innovation + 1
    return pool.innovation
end

-- Create a new pool of organisms
function newPool()
    local pool = {}
    pool.species = {}
    pool.generation = 0
    pool.innovation = outputCount
    pool.currentSpecies = 1
    pool.currentGenome = 1
    pool.currentFrame = 0
    pool.maxFitness = 0

    return pool
end

-- Create a new species
function newSpecies()
    local species = {}
    species.topFitness = 0
    species.staleness = 0
    species.genomes = {}
    species.averageFitness = 0

    return species
end

-- Create a new genome (practically a new organism/neural network)
function createGenome()
    local genome = {}
    genome.genes = {}
    genome.fitness = 0
    genome.adjustedFitness = 0
    genome.network = {}
    genome.maxneuron = 0
    genome.globalRank = 0
    genome.mutationRates = {}
    genome.mutationRates["connections"] = MutateConnectionsChance
    genome.mutationRates["link"] = LinkMutationChance
    genome.mutationRates["bias"] = BiasMutationChance
    genome.mutationRates["node"] = NodeMutationChance
    genome.mutationRates["enable"] = EnableMutationChance
    genome.mutationRates["disable"] = DisableMutationChance
    genome.mutationRates["step"] = StepSize

    return genome
end

-- Copy information from one genome to another
function copyGenome(genome)
    local genome2 = createGenome()
    for g = 1, #genome.genes do
        table.insert(genome2.genes, copyGene(genome.genes[g]))
    end
    genome2.maxneuron = genome.maxneuron
    genome2.mutationRates["connections"] = genome.mutationRates["connections"]
    genome2.mutationRates["link"] = genome.mutationRates["link"]
    genome2.mutationRates["bias"] = genome.mutationRates["bias"]
    genome2.mutationRates["node"] = genome.mutationRates["node"]
    genome2.mutationRates["enable"] = genome.mutationRates["enable"]
    genome2.mutationRates["disable"] = genome.mutationRates["disable"]

    return genome2
end

-- Create a blank/empty genome
function basicGenome()
    local genome = createGenome()
    local innovation = 1

    genome.maxneuron = inputCount
    mutate(genome)

    return genome
end

-- Create a blank/empty gene
function newGene()
    local gene = {}
    gene.inNode = 0
    gene.outNode = 0
    gene.weight = 0
    gene.expressed = true
    gene.innovation = 0

    return gene
end

-- Generate a copy of a gene
function copyGene(gene)
    local gene2 = newGene()
    gene2.inNode = gene.inNode
    gene2.outNode = gene.outNode
    gene2.weight = gene.weight
    gene2.expressed = gene.expressed
    gene2.innovation = gene.innovation

    return gene2
end

-- Create a new, blank neuron (input or output)
function newNeuron()
    local neuron = {}
    neuron.incoming = {}
    neuron.value = 0.0

    return neuron
end

-- Create a network of neurons and nodes
function generateNetwork(genome)
    local network = {}
    network.neurons = {}

    -- Create input neurons
    for i = 1, inputCount do
        network.neurons[i] = newNeuron()
    end

    -- Create output nodes
    for o = 1, outputCount do
        network.neurons[MaxNodes + o] = newNeuron()
    end

    table.sort(genome.genes, function(a, b)
        return (a.outNode < b.outNode)
    end)

    for i = 1, #genome.genes do
        local gene = genome.genes[i]
        if gene.expressed then
            if network.neurons[gene.outNode] == nil then
                network.neurons[gene.outNode] = newNeuron()
            end
            local neuron = network.neurons[gene.outNode]
            table.insert(neuron.incoming, gene)
            if network.neurons[gene.inNode] == nil then
                network.neurons[gene.inNode] = newNeuron()
            end
        end
    end

    genome.network = network
end

-- With a given network and inputs, return the proper outputs
function evaluateNetwork(network, inputs)
    table.insert(inputs, 1)

    -- Catch an error with the number of inputs
    if #inputs ~= inputCount then
        console.writeline("Incorrect number of neural network inputs.")
        return {}
    end

    for i = 1, inputCount do
        network.neurons[i].value = inputs[i]
    end

    for _, neuron in pairs(network.neurons) do
        local sum = 0
        for j = 1, #neuron.incoming do
            local incoming = neuron.incoming[j]
            local other = network.neurons[incoming.inNode]
            sum = sum + incoming.weight * other.value
        end

        if #neuron.incoming > 0 then
            neuron.value = sigmoid(sum)
        end
    end

    -- Set the outputs
    local outputs = {}
    for o = 1, outputCount do
        local button = "P1 " .. ButtonNames[o]
        if network.neurons[MaxNodes + o].value > 0 then
            outputs[button] = true
        else
            outputs[button] = false
        end
    end

    return outputs
end

-- Genetically crossover two genomes
function crossover(g1, g2)
    -- Make sure g1 is the higher fitness genome
    if g2.fitness > g1.fitness then
        g1, g2 = g2, g1
    end

    local child = createGenome()

    local innovations2 = {}
    for i = 1, #g2.genes do
        local gene = g2.genes[i]
        innovations2[gene.innovation] = gene
    end

    for i = 1, #g1.genes do
        local gene1 = g1.genes[i]
        local gene2 = innovations2[gene1.innovation]
        if gene2 ~= nil and math.random(2) == 1 and gene2.expressed then
            table.insert(child.genes, copyGene(gene2))
        else
            table.insert(child.genes, copyGene(gene1))
        end
    end

    child.maxneuron = math.max(g1.maxneuron, g2.maxneuron)

    for mutation, rate in pairs(g1.mutationRates) do
        child.mutationRates[mutation] = rate
    end

    return child
end

-- Pick a random neuron from a set of genes,
-- and say whether or not it should be an input node
function randomNeuron(genes, isInput)
    local neurons = {}
    if isInput then
        for i = 1, inputCount do
            neurons[i] = true
        end
    end
    for o = 1, outputCount do
        neurons[MaxNodes + o] = true
    end
    for i = 1, #genes do
        if (isInput) or genes[i].inNode > inputCount then
            neurons[genes[i].inNode] = true
        end
        if (isInput) or genes[i].outNode > inputCount then
            neurons[genes[i].outNode] = true
        end
    end

    local count = 0
    for _, _ in pairs(neurons) do
        count = count + 1
    end
    local n = math.random(1, count)

    for k, v in pairs(neurons) do
        n = n - 1
        if n == 0 then
            return k
        end
    end

    return 0
end


-- Point mutation for a given genome
function pointMutate(genome)
    local step = genome.mutationRates["step"]

    for i = 1, #genome.genes do
        local gene = genome.genes[i]
        if math.random() < PerturbChance then
            gene.weight = gene.weight + math.random() * step * 2 - step
        else
            gene.weight = math.random() * 4 - 2
        end
    end
end

-- This mutation takes two unconnected nodes, and adds a connection gene with a random weight between them.
function addConnectionMutation(genome, forceBias)
    local node1 = randomNeuron(genome.genes, true)
    local node2 = randomNeuron(genome.genes, false)
    local newLink = newGene()

    -- return if both nodes are input nodes
    if node2 <= inputCount and node1 <= inputCount then
        return
    end

    -- Swap output and input, then add the new link between the nodes
    if node1 <= inputCount then
        node1, node2 = node2, node1
    end
    newLink.inNode = node2
    newLink.outNode = node1
    if forceBias then
        newLink.inNode = inputCount
    end


    -- check to see if nodes are already connected
    for i = 1, #genome.genes do
        if genome.genes[i].inNode == newLink.inNode and genome.genes[i].outNode == newLink.outNode then
            return
        end
    end


    newLink.innovation = newInnovation()
    newLink.weight = math.random() * 4 - 2

    table.insert(genome.genes, newLink)
end

-- Node mutation for a given genome
function nodeMutate(genome)
    if #genome.genes == 0 then
        return
    end

    genome.maxneuron = genome.maxneuron + 1

    local gene = genome.genes[math.random(1, #genome.genes)]
    if not gene.expressed then
        return
    end
    gene.expressed = false

    local gene1 = copyGene(gene)
    gene1.outNode = genome.maxneuron
    gene1.weight = 1.0
    gene1.innovation = newInnovation()
    gene1.expressed = true
    table.insert(genome.genes, gene1)

    local gene2 = copyGene(gene)
    gene2.inNode = genome.maxneuron
    gene2.innovation = newInnovation()
    gene2.expressed = true
    table.insert(genome.genes, gene2)
end

-- Swap the state of enabled or disabled mutation
function enableDisableMutate(genome, enable)
    local candidates = {}
    for _, gene in pairs(genome.genes) do
        if gene.expressed == not enable then
            table.insert(candidates, gene)
        end
    end

    if #candidates == 0 then
        return
    end

    local gene = candidates[math.random(1, #candidates)]
    gene.expressed = not gene.expressed
end


-- Handle all kinds of mutation for a given genome
function mutate(genome)
    for mutation, rate in pairs(genome.mutationRates) do
        if math.random(1, 2) == 1 then
            genome.mutationRates[mutation] = 0.95 * rate
        else
            genome.mutationRates[mutation] = 1.05263 * rate
        end
    end

    if math.random() < genome.mutationRates["connections"] then
        pointMutate(genome)
    end

    local p = genome.mutationRates["link"]
    while p > 0 do
        if math.random() < p then
            addConnectionMutation(genome, false)
        end
        p = p - 1
    end

    p = genome.mutationRates["bias"]
    while p > 0 do
        if math.random() < p then
            addConnectionMutation(genome, true)
        end
        p = p - 1
    end

    p = genome.mutationRates["node"]
    while p > 0 do
        if math.random() < p then
            nodeMutate(genome)
        end
        p = p - 1
    end

    p = genome.mutationRates["enable"]
    while p > 0 do
        if math.random() < p then
            enableDisableMutate(genome, true)
        end
        p = p - 1
    end

    p = genome.mutationRates["disable"]
    while p > 0 do
        if math.random() < p then
            enableDisableMutate(genome, false)
        end
        p = p - 1
    end
end

-- Check for now many disjoint genes there are in a pair
function disjoint(genes1, genes2)
    local i1 = {}
    for i = 1, #genes1 do
        local gene = genes1[i]
        i1[gene.innovation] = true
    end

    local i2 = {}
    for i = 1, #genes2 do
        local gene = genes2[i]
        i2[gene.innovation] = true
    end

    local disjointGenes = 0
    for i = 1, #genes1 do
        local gene = genes1[i]
        if not i2[gene.innovation] then
            disjointGenes = disjointGenes + 1
        end
    end

    for i = 1, #genes2 do
        local gene = genes2[i]
        if not i1[gene.innovation] then
            disjointGenes = disjointGenes + 1
        end
    end

    local n = math.max(#genes1, #genes2)

    return disjointGenes / n
end

-- Find the average weights between two genes
function weights(genes1, genes2)
    local i2 = {}
    for i = 1, #genes2 do
        local gene = genes2[i]
        i2[gene.innovation] = gene
    end

    local sum = 0
    local coincident = 0
    for i = 1, #genes1 do
        local gene = genes1[i]
        if i2[gene.innovation] ~= nil then
            local gene2 = i2[gene.innovation]
            sum = sum + math.abs(gene.weight - gene2.weight)
            coincident = coincident + 1
        end
    end

    return sum / coincident
end

-- See if two genomes can fit into the same species,
-- based on how similar their neural networks are
function sameSpecies(genome1, genome2)
    local dd = DeltaDisjoint * disjoint(genome1.genes, genome2.genes)
    local dw = DeltaWeights * weights(genome1.genes, genome2.genes)
    return dd + dw < DeltaThreshold
end

-- Rank all of the genomes by fitness
function rankGlobally()
    local global = {}
    for s = 1, #pool.species do
        local species = pool.species[s]
        for g = 1, #species.genomes do
            table.insert(global, species.genomes[g])
        end
    end
    table.sort(global, function(a, b)
        return (a.fitness < b.fitness)
    end)

    for g = 1, #global do
        global[g].globalRank = g
    end
end

-- Return the average fitness in a species
function calculateAverageFitness(species)
    local total = 0

    for g = 1, #species.genomes do
        local genome = species.genomes[g]
        total = total + genome.globalRank
    end

    species.averageFitness = total / #species.genomes
end

-- Total the average fitnesses in a species
function totalAverageFitness()
    local total = 0
    for s = 1, #pool.species do
        local species = pool.species[s]
        total = total + species.averageFitness
    end

    return total
end

-- Remove the weaker genomes in each species
function cullSpecies(cutToOne)
    for s = 1, #pool.species do
        local species = pool.species[s]

        table.sort(species.genomes, function(a, b)
            return (a.fitness > b.fitness)
        end)

        local remaining = math.ceil(#species.genomes / 2)
        if cutToOne then
            remaining = 1
        end
        while #species.genomes > remaining do
            table.remove(species.genomes)
        end
    end
end

-- Crossover two genomes from a species, and then mutate
function breedChild(species)
    local child = {}
    if math.random() < CrossoverChance then
        g1 = species.genomes[math.random(1, #species.genomes)]
        g2 = species.genomes[math.random(1, #species.genomes)]
        child = crossover(g1, g2)
    else
        g = species.genomes[math.random(1, #species.genomes)]
        child = copyGenome(g)
    end

    mutate(child)

    return child
end

-- Remove the species that have stayed around for too long
function removeStaleSpecies()
    local survived = {}

    for s = 1, #pool.species do
        local species = pool.species[s]

        table.sort(species.genomes, function(a, b)
            return (a.fitness > b.fitness)
        end)

        if species.genomes[1].fitness > species.topFitness then
            species.topFitness = species.genomes[1].fitness
            species.staleness = 0
        else
            species.staleness = species.staleness + 1
        end
        if species.staleness < StaleSpecies or species.topFitness >= pool.maxFitness then
            table.insert(survived, species)
        end
    end

    pool.species = survived
end

-- Remove all of the weak species
function removeWeakSpecies()
    local survived = {}

    local sum = totalAverageFitness()
    for s = 1, #pool.species do
        local species = pool.species[s]
        breed = math.floor(species.averageFitness / sum * Population)
        if breed >= 1 then
            table.insert(survived, species)
        end
    end

    pool.species = survived
end

-- Find a species for the offspring, or create a new one
function addToSpecies(child)
    local foundSpecies = false
    for s = 1, #pool.species do
        local species = pool.species[s]
        if not foundSpecies and sameSpecies(child, species.genomes[1]) then
            table.insert(species.genomes, child)
            foundSpecies = true
        end
    end

    -- Handle the new child if they don't have a species
    if not foundSpecies then
        local childSpecies = newSpecies()
        table.insert(childSpecies.genomes, child)
        table.insert(pool.species, childSpecies)
    end
end

-- Handle everything to create a new generation
function newGeneration()
    cullSpecies(false) -- Remove the bottom half of each species
    rankGlobally() -- Rank the genomes
    removeStaleSpecies() -- Remove the necessary species
    rankGlobally() -- Re-rank them
    -- Calculate each species' average fitness
    for s = 1, #pool.species do
        local species = pool.species[s]
        calculateAverageFitness(species)
    end
    removeWeakSpecies() -- Remove the weak species after calculation
    local sum = totalAverageFitness()
    local children = {}
    for s = 1, #pool.species do
        local species = pool.species[s]
        breed = math.floor(species.averageFitness / sum * Population) - 1 -- Determine how many children to make
        for i = 1, breed do
            table.insert(children, breedChild(species))
        end
    end
    cullSpecies(true) -- Remove all but the top member of each species
    -- Keep creating children until the pool is full
    while #children + #pool.species < Population do
        local species = pool.species[math.random(1, #pool.species)]
        table.insert(children, breedChild(species))
    end
    -- Add the new children to a species
    for c = 1, #children do
        local child = children[c]
        addToSpecies(child)
    end

    -- Increment generation
    pool.generation = pool.generation + 1

--    writeFile("backup." .. pool.generation .. "." .. forms.gettext(saveLoadFile))
end

-- Create a new pool
function initializePool()
    pool = newPool()

    for i = 1, Population do
        basic = basicGenome()
        addToSpecies(basic)
    end

    initializeRun()
end

-- Make sure that the input nodes (and controller) are currently empty
function clearJoypad()
    controller = {}
    for b = 1, #ButtonNames do
        controller["P1 " .. ButtonNames[b]] = false
    end
    joypad.set(controller)
end

-- Start the game with a new generation
function initializeRun()
    savestate.load(Filename);
    bestScore = 0
    pool.currentFrame = 0
    timeout = 60
    clearJoypad()

    local species = pool.species[pool.currentSpecies]
    local genome = species.genomes[pool.currentGenome]
    generateNetwork(genome)
    evaluateCurrent()
end

-- Evaluate the current input
function evaluateCurrent()
    local species = pool.species[pool.currentSpecies]
    local genome = species.genomes[pool.currentGenome]

    -- Get the inputs and handle them accordingly
    inputs = genInputs()
    controller = evaluateNetwork(genome.network, inputs)

    -- Account for both left and right being pressed
    if controller["P1 Left"] and controller["P1 Right"] then
        controller["P1 Left"] = false
        controller["P1 Right"] = false
    end

    -- Same thing for up and down
    if controller["P1 Up"] and controller["P1 Down"] then
        controller["P1 Up"] = false
        controller["P1 Down"] = false
    end

    joypad.set(controller)
end

-- If there is no pool, then create a new one
if pool == nil then
    initializePool()
end

-- Iterate to the next genome (or species, or generation)
function nextGenome()
    pool.currentGenome = pool.currentGenome + 1
    if pool.currentGenome > #pool.species[pool.currentSpecies].genomes then
        pool.currentGenome = 1
        pool.currentSpecies = pool.currentSpecies + 1
        if pool.currentSpecies > #pool.species then
            newGeneration()
            pool.currentSpecies = 1
        end
    end
end

-- See if the fitness has already been scored
function fitnessAlreadyMeasured()
    local species = pool.species[pool.currentSpecies]
    local genome = species.genomes[pool.currentGenome]

    return genome.fitness ~= 0
end




-- Find the absolute strongest organism
function playTop()

    -- Iterate through the organisms
    local maxfitness = 0
    local maxs, maxg
    for s, species in pairs(pool.species) do
        for g, genome in pairs(species.genomes) do
            -- Take action if they're stronger than the last
            if genome.fitness > maxfitness then
                maxfitness = genome.fitness
                maxs = s
                maxg = g
            end
        end
    end

    -- Perform actions with the strongest
    pool.currentSpecies = maxs
    pool.currentGenome = maxg
    pool.maxFitness = maxfitness
--    forms.settext(maxFitnessLabel, "Max Fitness: " .. math.floor(pool.maxFitness))
    initializeRun()
    pool.currentFrame = pool.currentFrame + 1
    return
end

-- Exit the game
--function onExit()
--    forms.destroy(form)
--end

-- Write to file, handle exit sequence, and draw information on the screen
--writeFile("temp.pool")

event.onexit(onExit)

--form = forms.newform(200, 260, "Fitness")
--maxFitnessLabel = forms.label(form, "Max Fitness: " .. math.floor(pool.maxFitness), 5, 8)
--playTopButton = forms.button(form, "Play Top", playTop, 5, 170)

-- Constantly update the game
while true do

    -- Update the species and genome
    local species = pool.species[pool.currentSpecies]
    local genome = species.genomes[pool.currentGenome]

    -- evaluate frames 10 times per second.  NES runs at 60fps
    if pool.currentFrame % 6 == 0 then
        evaluateCurrent()
    end

    joypad.set(controller)

    updateCharacterPositions()
    if 192 - pelletCount > bestScore then
        bestScore = 192 - pelletCount
        timeout = 20 + 10 * (192 - pelletCount)
        --        timeout = 60
    end

    timeout = timeout - 1

    gui.drawText(175, 90, "gen:" .. pool.generation)
    gui.drawText(175, 105, "species:" .. pool.currentSpecies)
    if(pool.currentGenome > 1) then
        gui.drawText(175, 120, "genome:" .. pool.currentGenome -1)
    else
        gui.drawText(175, 120, "genome:" .. pool.currentGenome)
    end
    gui.drawText(175, 150, "top:" .. pool.maxFitness)
    gui.drawText(175, 165, "Q - replay")

    local timeoutBonus = pool.currentFrame / 2
    if timeout + timeoutBonus <= 0 or lifeCount < 3 then

        local fitness = bestScore * 10

        if fitness == 0 then
            fitness = -1
        end
        genome.fitness = fitness

        -- If it's a new record, then record it and write it to the file
        if fitness > pool.maxFitness then
            pool.maxFitness = fitness
--            forms.settext(maxFitnessLabel, "Max Fitness: " .. math.floor(pool.maxFitness))
        end

--        console.writeline("Gen " .. pool.generation .. " species " .. pool.currentSpecies .. " genome " .. pool.currentGenome .. " fitness: " .. fitness)
        pool.currentSpecies = 1
        pool.currentGenome = 1


        -- Continue iterating through genomes that have already been scored for fitness
        while fitnessAlreadyMeasured() do
            nextGenome()
        end
        initializeRun()
    end

    -- Total the number of measured genomes
    local measured = 0
    local total = 0
    for _, species in pairs(pool.species) do
        for _, genome in pairs(species.genomes) do
            total = total + 1
            if genome.fitness ~= 0 then
                measured = measured + 1
            end
        end
    end


    local keys = input.get()
    if(keys['Q'] == true) then
        playTop()
    end

    pool.currentFrame = pool.currentFrame + 1
    emu.frameadvance();
end