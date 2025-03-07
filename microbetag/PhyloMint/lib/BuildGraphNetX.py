#!/usr/bin/env python3
from libsbml import readSBML
import networkx as nx


def buildDG(sbml):
    '''
    Usage: reads SBML file, parses reaction and product list
    Returns: networkx directed graph
    '''

    # initate empty directed graph
    DG = nx.DiGraph()

    document = readSBML(sbml)
    model = document.getModel()

    # # get reactions
    # rxns = (model.getListOfReactions())

    # get exchange reactions
    cellular_reactions = [r for r in model.getListOfReactions() if not r.getId().startswith("R_EX_")]
    exchange_reactions = [r for r in model.getListOfReactions() if r.getId().startswith("R_EX_")]

    # get reaction and product and construct directed graph
    # for rxn in rxns:
    for rxn in cellular_reactions:

        react = [i.getSpecies() for i in rxn.getListOfReactants()]
        prod = [j.getSpecies() for j in rxn.getListOfProducts()]
        for r in react:
            for p in prod:
                DG.add_edge(r,p)

        # NOTE (2025-03-07): In case of reversible reactions, we consider that too
        if rxn.reversible:
            react = [i.getSpecies() for i in rxn.getListOfProducts()]
            prod = [j.getSpecies() for j in rxn.getListOfReactants()]
            for r in react:
                for p in prod:
                    DG.add_edge(r,p)
    return DG

def getSeedSet(DG, maxComponentSize = 5):
    '''
    Usage: takes input networkX directed graph
    Returns: SeedSet dictionary{seedset:confidence score}
    Implementation follows literature description,
    Improves upon NetCooperate module implementation which erroneously discards certian cases of SCCs (where a smaller potential SCC lies within a larger SCC)
    '''

    # get SCC
    SCC = nx.strongly_connected_components(DG)

    SeedSetConfidence = dict()

    for cc in SCC:
        # convert set to list
        cc_temp = list(cc)

        # filter out CC larger than threshold
        if len(cc_temp) > maxComponentSize:
            continue

        # check single element SCC
        elif len(cc_temp) == 1:
            if DG.in_degree(cc_temp[0]) == 0:
                SeedSetConfidence[cc_temp[0]] = 1.0

        # check 2 to max threshold SCC
        else:
            #check if no out nodes
            for node in cc_temp:
                # check every edge of SCC
                for edge in DG.in_edges(node):
                    # if SCC is not self contained, then it is not considered seed set
                    if edge[0] not in cc_temp:
                        cc_temp = []
            for node in cc_temp:
                SeedSetConfidence[node] = 1/len(cc_temp)

    SeedSet = SeedSetConfidence.keys()

    nonSeedSet = list(set(DG.nodes()) - set(SeedSet))

    return(SeedSetConfidence, SeedSet, nonSeedSet)
