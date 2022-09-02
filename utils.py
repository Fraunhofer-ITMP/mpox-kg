# -*- coding: utf-8 -*-

"""Utils files with all functions relevant to generation of KG."""

import logging
import pickle
from collections import defaultdict

import pandas as pd
import pubchempy as pcp
import requests
from chembl_webresource_client.http_errors import HttpBadRequest, HttpApplicationError
from chembl_webresource_client.new_client import new_client
from pybel import BELGraph
from pybel.dsl import Protein, Abundance, Pathology, BiologicalProcess
from tqdm import tqdm
import time

logger = logging.getLogger("__name__")


def RetMech(chemblIds) -> dict:
    """Function to retrieve mechanism of actions and target proteins from ChEMBL

    :param chemblIds:
    :return:
    """
    getMech = new_client.mechanism

    mechList = []
    for chemblid in tqdm(chemblIds, desc='Retrieving mechanisms from ChEMBL'):
        mechs = getMech.filter(
            molecule_chembl_id=chemblid
        ).only(['mechanism_of_action', 'target_chembl_id'])
        mechList.append(list(mechs))

    named_mechList = dict(zip(chemblIds, mechList))
    named_mechList = {
        k: v
        for k, v in named_mechList.items()
        if v
    }
    return named_mechList


def RetDrugInd(chemblIDs) -> dict:
    """Function to retrieve associated diseases from ChEMBL

    :param chemblIDs:
    :return:
    """
    getDrugInd = new_client.drug_indication

    drugIndList = []
    for chemblid in tqdm(chemblIDs, desc='Retrieving diseases from ChEMBL'):
        drugInd = getDrugInd.filter(
            molecule_chembl_id=chemblid
        ).only('mesh_heading')
        drugIndList.append(list(drugInd))

    named_drugIndList = dict(zip(chemblIDs, drugIndList))
    named_drugIndList = {
        k: v
        for k, v in named_drugIndList.items()
        if v
    }
    return named_drugIndList

def RetAct(chemblIds) -> dict:
    """Function to retrieve associated assays from ChEMBL

    :param chemblIds:
    :return:
    """
    GetAct = new_client.activity
    ActList = []
    filtered_list=['assay_chembl_id','assay_type','pchembl_value','target_chembl_id',
                   'target_organism','bao_label','target_type']

#     filtered_list = [
#         'pchembl_value',
#         'target_chembl_id',
#         'target_type',
#         'bao_label'
#     ]

    for chembl in tqdm(chemblIds, desc='Retrieving bioassays from ChEMBL'):
    #for i in range(len(chemblIds)):
        acts = GetAct.filter(
            molecule_chembl_id=chembl,
            pchembl_value__isnull=False,
            assay_type_iregex='(B|F)',
            target_organism='Homo sapiens'
        ).only(filtered_list)
        
        #print(chemblIds[i])
        data = []

        for d in acts:

            if float(d.get('pchembl_value')) < 6:
                continue
            
#             try:
#                 if d['target_type'] in ('CELL-LINE', 'UNCHECKED'):
#                     continue
#             except KeyError:
#                 continue
            if (d.get('bao_label') != 'single protein format'):
                continue

            #uprot_id = d['target_components'][0]['accession']
            #print(uprot_id)
            #print(d)
            data.append(d)

        # acts = acts[:5]
        ActList.append(list(data))

    named_ActList = dict(zip(chemblIds, ActList))
    named_ActList = {
        k: v
        for k, v in named_ActList.items()
        if v
    }
    return named_ActList

# def RetAct(chemblIds) -> dict:
    # """Function to retrieve associated assays from ChEMBL

    # :param chemblIds:
    # :return:
    # """
    # GetAct = new_client.activity
    # ActList = []

    # filtered_list = [
        # 'pchembl_value',
        # 'target_chembl_id',
        # 'target_type'
    # ]

    # for chembl in tqdm(chemblIds, desc='Retrieving bioassays from ChEMBL'):
        # acts = GetAct.filter(
            # molecule_chembl_id=chembl,
            # pchembl_value__isnull=False,
            # assay_type_iregex='(B|F)',
            # target_organism='Homo sapiens'
        # ).only(filtered_list)

        # data = []

        # for d in acts:

            # if float(d.get('pchembl_value')) < 6:
                # continue
            # try:

                # if d['target_type'] in ('CELL-LINE', 'UNCHECKED'):
                    # continue
            # except IndexError:
                # continue

            # uprot_id = d['target_components'][0]['accession']

            # data.append(uprot_id)

        # # acts = acts[:5]
        # ActList.append(list(data))

    # named_ActList = dict(zip(chemblIds, ActList))
    # named_ActList = {
        # k: v
        # for k, v in named_ActList.items()
        # if v
    # }
    # return named_ActList


def Ret_chembl_protein(sourceList) -> list:
    """Method to retrieve ChEMBL ids which are proteins/targets

    :param sourceList:
    :return:
    """
    protein_List = []
    for item in sourceList:
        for j in range(len(sourceList[item])):
            protein_List.append(sourceList[item][j]['target_chembl_id'])

    protein_List = set(protein_List)
    protein_List = list(filter(None, protein_List))
    return protein_List


def chembl2uniprot(chemblIDs) -> dict:
    """Method to convert ChEMBL id to UNIPROT and get associated REACTOME pathways

    :param chemblIDs:
    :return:
    """
    getTarget = new_client.target
    chem2Gene2path = []
    chemHasNoPath = set()
    chemNotprotein = set()

    chem2path = defaultdict(list)

    # Loop to ensure it is a protein
    for chemblid in tqdm(chemblIDs, desc='Filtering UniProt proteins from ChEMBL'):
        chem = getTarget.filter(
            chembl_id=chemblid
        ).only('target_components')

        try:
            uprot_id = chem[0]['target_components'][0]['accession']

            if not uprot_id:
                chemHasNoPath.add(chemblid)

        except IndexError:
            chemHasNoPath.add(chemblid)

    logger.info(f'No UniProt information available for {len(chemHasNoPath)} proteins.')

    chemblIDs_filtered = [
        item
        for item in chemblIDs
        if item not in chemHasNoPath
    ]

    # Get gene symbol from ChEMBL and filtering the list for human proteins only
    for chemblid in tqdm(chemblIDs_filtered, desc='Filtering human proteins from ChEMBL'):

        chem = getTarget.filter(chembl_id=chemblid).only('target_components')
        getGene = chem[0]['target_components'][0]['target_component_synonyms']
        try:
            getGene = [item for item in getGene if item["syn_type"] == "GENE_SYMBOL"][0]

            if not getGene:
                chemNotprotein.add(chemblid)

        except IndexError:
            chemNotprotein.add(chemblid)

    chemblIDs_filtered = [
        item
        for item in chemblIDs_filtered
        if item not in chemNotprotein
    ]

    # Extracting data for valid proteins only
    for chemblid in tqdm(chemblIDs_filtered, desc='Populating ChEMBL data for human proteins'):
        chem = getTarget.filter(
            chembl_id=chemblid
        ).only('target_components')

        # UniProt data
        uprot_id = chem[0]['target_components'][0]['accession']

        # Gene symbol
        getGene = chem[0]['target_components'][0]['target_component_synonyms']
        getGene = [item for item in getGene if item["syn_type"] == "GENE_SYMBOL"][0]

        # Pathway data
        chem2path = [item for item in chem[0]['target_components'][0]['target_component_xrefs'] if
                     item["xref_src_db"] == "Reactome"]

        uprot = {'accession': uprot_id}
        chem2path.append(uprot)
        chem2path.append(getGene)
        chem2Gene2path.append(chem2path)

    named_chem2Gene2path = dict(zip(chemblIDs_filtered, chem2Gene2path))
    named_chem2Gene2path = {
        k: v
        for k, v in named_chem2Gene2path.items()
        if v
    }
    return named_chem2Gene2path


def chembl2gene2path(
    chem2geneList,
    ActList
):
    """Method for updating chembl protein nodes with gene symbol.

    :param chem2geneList:
    :param ActList:
    :return:
    """
    for item in chem2geneList:
        sizeOfitem = len(chem2geneList[item])
        gene = chem2geneList[item][sizeOfitem - 1]['component_synonym']
        for jtem in ActList:
            for i in range(len(ActList[jtem])):
                if item == ActList.get(jtem)[i]['target_chembl_id']:
                    newkey = {'Protein': gene}
                    ActList[jtem][i].update(newkey)

    return ActList


def Ret_uprotid(chembl2uprot) -> list:
    """Method to get UniProt ids from dict of chembl2uniprot

    :param chembl2uprot:
    :return:
    """
    chemblProt = []
    for item in chembl2uprot:
        for j in range(len(chembl2uprot[item])):
            if 'accession' in chembl2uprot[item][j]:
                chemblProt.append(
                    chembl2uprot[item][j]['accession']
                )
    return chemblProt


def ExtractFromUniProt(uniprot_id) -> dict:
    """Uniprot parser to retrieve information about OMIM disease, reactome pathway, biological process,
     and molecular functions.

    :param uniprot_id:
    :return:
    """
    Uniprot_Dict = []

    mapped_uprot = []

    for id in uniprot_id:

        # Retrieve data for id in text format if found in uniprot
        ret_uprot = requests.get(
            'https://www.uniprot.org/uniprot/' + id + '.txt'
        ).text.split('\n')

        if ret_uprot == ['']:
            continue

        id_copy = id
        mapped_uprot.append(id_copy)
        i = 0
        j = 0
        k = 0
        id = {}
        id['Disease'] = {}
        id['Reactome'] = {}
        id['Function'] = {}
        id['BioProcess'] = {}
        id['Gene'] = {}

        # parse each line looking for info about disease, pathway, funcn, bp and so on
        for line in ret_uprot:
            # parse lines with disease and extract disease names and omim ids
            if '-!- DISEASE:' in line:
                if ('[MIM:' in line):
                    dis = line.split(':')
                    id['Disease'].update({dis[1][1:-5]: dis[2][:-1]})

            # extract reactome ids and names
            if 'Reactome;' in line:
                ract = line.split(';')
                id['Reactome'].update({ract[2][1:-2]: ract[1][1:]})

            # look for functions
            if ' F:' in line:
                if j < 5:
                    fn = line.split(';')
                    id['Function'].update({fn[2][3:]: fn[1][1:]})
                    j += 1

            # look for biological processes
            if ' P:' in line:
                if i < 5:
                    bp = line.split(';')
                    # bp returns list with GO ids and names
                    id['BioProcess'].update({bp[2][3:]: bp[1][1:]})
                    i += 1

            if 'GN   Name' in line:
                if k == 0:
                    gene = line.split('=')
                    gene = gene[1].split(' ')
                    if ';' in gene[0]:
                        gene = gene[0].split(';')
                        gene = {'Gene': gene[0]}
                    else:
                        gene = {'Gene': gene[0]}
                    id.update(gene)
                    k += 1

        Uniprot_Dict.append(id)

    Uniprot_Dict = dict(zip(mapped_uprot, Uniprot_Dict))

    return Uniprot_Dict


def chem2moa_rel(
    named_mechList,
    org,
    graph: BELGraph
) -> BELGraph:
    """Method to create the monkeypox graph

    :param named_mechList:
    :param org:
    :param graph: BEL graph of Monkeypox
    :return:
    """
    for chembl_name, chembl_entries in tqdm(named_mechList.items(), desc='Populating Chemical-MoA edges'):
        for info in chembl_entries:
            graph.add_association(
                Abundance(namespace='ChEMBL', name=chembl_name),
                BiologicalProcess(namespace='MOA', name=info['mechanism_of_action']),
                citation='ChEMBL database',
                evidence='ChEMBL query'
            )

            if not info['target_chembl_id']:
                continue

            if 'Protein' in info:
                graph.add_association(
                    Abundance(namespace='ChEMBL', name=chembl_name),
                    Protein(namespace=org, name=info['Protein']),
                    citation='ChEMBL database',
                    evidence='ChEMBL query'
                )
            else:
                graph.add_association(
                    Abundance(namespace='ChEMBL', name=chembl_name),
                    Protein(namespace=org, name=info['target_chembl_id']),
                    citation='ChEMBL database',
                    evidence='ChEMBL query'
                )

    return graph


def chem2dis_rel(
    named_drugIndList,
    graph: BELGraph
) -> BELGraph:
    """Method to add drug indication edges to the KG.

    :param named_drugIndList:
    :param graph:
    :return:
    """
    for chembl_id, drug_entries in tqdm(named_drugIndList.items(), desc='Populating Drug-Indication edges'):
        for drug_data in drug_entries:
            graph.add_association(
                Abundance(namespace='ChEMBL', name=chembl_id),
                Pathology(namespace='Disease', name=drug_data['mesh_heading']),
                citation='ChEMBL database',
                evidence='ChEMBL query'
            )
    return graph


def chem2act_rel(
    named_ActList,
    org,
    graph: BELGraph
) -> BELGraph:
    """Method to add bioassay edges to the KG.

    :param named_ActList:
    :param org:
    :param graph:
    :return:
    """
    for chemical, chem_entries in tqdm(named_ActList.items(), desc='Adding bioassay edges to BEL'):
        for chem_data in chem_entries:
            if chem_data['target_chembl_id']:
                if 'Protein' in chem_data:
                    graph.add_association(
                        Abundance(namespace='ChEMBLAssay', name=chem_data['assay_chembl_id']),
                        Protein(namespace=org, name=chem_data['Protein']),
                        citation='ChEMBL database',
                        evidence='ChEMBL query'
                    )
                else:
                    graph.add_association(
                        Abundance(namespace='ChEMBLAssay', name=chem_data['assay_chembl_id']),
                        Protein(namespace=org, name=chem_data['target_chembl_id']),
                        citation='ChEMBL database',
                        evidence='ChEMBL query'
                    )

            graph.add_association(
                Abundance(namespace='ChEMBL', name=chemical),
                Abundance(namespace='ChEMBLAssay', name=chem_data['assay_chembl_id']),
                citation='ChEMBL database',
                evidence='ChEMBL query',
                annotation={
                    'assayType': chem_data['assay_type'],
                    'pChEMBL': chem_data['pchembl_value']
                }
            )

    return graph


def gene2path_rel(
    named_chem2geneList,
    org,
    graph
) -> BELGraph:
    """Method to add protein and reactome data to KG

    :param named_chem2geneList:
    :param org:
    :param graph:
    :return:
    """
    for item in named_chem2geneList:
        itemLen = len(named_chem2geneList[item]) - 1
        for j in range(itemLen - 1):
            graph.add_association(
                Protein(namespace=org, name=named_chem2geneList[item][itemLen]['component_synonym']),
                BiologicalProcess(namespace='Reactome', name=named_chem2geneList[item][j]['xref_name']),
                citation='ChEMBL database',
                evidence='ChEMBL query',
                annotation={
                    'Reactome': named_chem2geneList[item][j]['xref_id']
                }
            )

    return graph


def uniprot_rel(
    named_uprotList,
    org,
    graph
) -> BELGraph:
    """Method to add UniProt related edges

    :param named_uprotList:
    :param org:
    :param graph:
    :return:
    """
    for item in named_uprotList:
        fun = list(named_uprotList[item]['Function'].keys())
        bp = list(named_uprotList[item]['BioProcess'].keys())
        for f in fun:
            if str(named_uprotList[item]['Gene']) != 'nan' and not isinstance(named_uprotList[item]['Gene'], dict):
                graph.add_association(
                    Protein(namespace=org, name=named_uprotList[item]['Gene']),
                    BiologicalProcess(namespace='GOMF', name=f),
                    citation='UniProt database',
                    evidence='UniProt query'
                )
            else:
                graph.add_association(
                    Protein(namespace=org, name=item),
                    BiologicalProcess(namespace='GOMF', name=f),
                    citation='UniProt database',
                    evidence='UniProt query'
                )

        for b in bp:
            if str(named_uprotList[item]['Gene']) != 'nan' and not isinstance(named_uprotList[item]['Gene'], dict):
                graph.add_association(
                    Protein(namespace=org, name=named_uprotList[item]['Gene']),
                    BiologicalProcess(namespace='GOBP', name=b),
                    citation='UniProt database',
                    evidence='UniProt query'
                )
            else:
                graph.add_association(
                    Protein(namespace=org, name=item),
                    BiologicalProcess(namespace='GOBP', name=b),
                    citation='UniProt database',
                    evidence='UniProt query'
                )

    return graph


def _get_target_data(
    protein_list: list,
    organism: str
) -> pd.DataFrame:
    """Get chemical for target data from ChEMBL.

    :param protein_list:
    :param organism:
    :return:
    """
    df_data = []

    target = new_client.target
    activity = new_client.activity

    for protein in tqdm(protein_list, desc='Retrieving chemicals for proteins'):
        if pd.isna(protein):
            continue
        try:
            prot_data = [target.search(protein)[0]]

            # Search for protein with same synonym
            if prot_data == [None]:
                prot_data = target.filter(
                    target_synonym__icontains=protein, target_organism__istartswith=organism
                ).only(['target_chembl_id', 'target_pref_name', 'molecule_chembl_id', 'molecule_pref_name'])
        except HttpBadRequest:
            print(f'Cannot search for {protein} due to chembl error')
            continue

        # No results found
        if not prot_data:
            continue

        for prot in prot_data:
            # Absence of chembl id
            if not prot['target_chembl_id']:
                continue

            prot_activity_data = activity.filter(
                target_chembl_id=prot['target_chembl_id'],
                assay_type_iregex='(B|F)',
            ).only([
                'pchembl_value', 'molecule_chembl_id', 'activity_id', 'target_pref_name', 'molecule_pref_name'
            ])

            if len(prot_activity_data) < 1:
                continue

            for i in prot_activity_data:
                tmp = {}

                try:
                    if i['pchembl_value'] is None:
                        continue
                except HttpApplicationError:
                    continue

                pchembl_val = i['pchembl_value']

                if float(pchembl_val) < 6:
                    tmp['activity'] = 'inhibitor'
                else:
                    tmp['activity'] = 'activator'

                tmp.update({
                    'protein_symbol': protein,
                    'protein_name': i['target_pref_name'],
                    'aid': str(i['activity_id']),
                    'chembl_id': i['molecule_chembl_id'],
                    'compound_name': i['molecule_pref_name'].capitalize() if i['molecule_pref_name'] else ''
                })
                df_data.append(tmp)
            
            #time.sleep(1)
            #print('1 seconds of break')

    # Merge duplicated protein-chemical entries into one
    df = pd.DataFrame()

    for idx, row in tqdm(enumerate(df_data), total=len(df_data), desc='Preparing data'):
        if idx == 0:
            df = df.append(row, ignore_index=True)
        else:
            _in_df = df.loc[
                (df['protein_symbol'] == row['protein_symbol']) & (df['chembl_id'] == row['chembl_id'])
                ]

            if _in_df.empty:
                df = df.append(row, ignore_index=True)
            else:
                row_index = _in_df.index

                # Check existing citations
                existing_assays = set(df.loc[row_index, 'aid'].values[0].split(' | '))
                old_count = len(existing_assays)
                existing_assays.add(row['aid'])
                new_count = len(existing_assays)

                # Check if new citation added, if yes - add respective data
                if old_count < new_count:
                    df.loc[row_index, 'aid'] = ' | '.join(existing_assays)
    df = df[['activity', 'protein_symbol', 'protein_name', 'aid', 'chembl_id', 'compound_name']]
    return df


def target_list_to_chemical(
    proteins: list,
    organism: str = 'Homo sapiens',
) -> pd.DataFrame:
    """Extract chemical information on list of targets
    Usage:
    >> target_list_to_chemical(proteins=['RIPK'])
    """

    df = _get_target_data(protein_list=proteins, organism=organism)
    return df

def getNodeList(nodeName,itmpGraph):
    node_list = []
    for node in itmpGraph.nodes():
        if isinstance(node,pybel.dsl.Abundance):
            if node.namespace == nodeName:
                node_list.append(node.name)
    return(node_list)

def chembl2rxn_rel(
    chemblid_list,
    graph: BELGraph
) -> BELGraph:
    """

    :param chemblid_list:
    :param graph:
    :return:
    """
    infile = open('data/opentargets/drugReactions.pkl', 'rb')
    rxn_df = pickle.load(infile)
    infile.close()

    chembl_id_rxn = rxn_df[rxn_df['chembl_id'].isin(chemblid_list)]
    chembl_id_rxn = chembl_id_rxn.reset_index(drop=True)
    for i in range(len(chembl_id_rxn)):
        graph.add_association(
            Abundance(namespace='ChEMBL', name=chembl_id_rxn['chembl_id'][i]),
            Pathology(namespace='SideEffect', name=chembl_id_rxn['event'][i]),  # TODO: Fix namespace
            citation="OpenTargets Platform",
            evidence='DrugReactions'
        )

    return graph


def cid2chembl(cidList) -> list:
    """Method to convert Pubchem CIDs to ChEMBL ids

    :param cidList:
    :return:
    """

    cid2chembl_list = []

    for id in tqdm(cidList, desc='Converting PubChem ids to ChEMBL ids'):
        c = pcp.Compound.from_cid(id)

        for synonym in c.synonyms:
            if synonym.startswith('CHEMBL'):
                cid2chembl_list.append(synonym)

    return cid2chembl_list
    
# def chembl2rxn_rel(itmpGraph):
    
    # infile = open('data/normalized_data/drugReactions.pkl','rb')
    # rxn_df = pickle.load(infile)
    # infile.close()
    
    # chembl_id = []
    # for node in itmpGraph.nodes():
        # if isinstance(node,pybel.dsl.Abundance):
            # if node.namespace == 'ChEMBL':
                # chembl_id.append(node.name)
            
    # chembl_id_rxn = rxn_df[rxn_df['chembl_id'].isin(chembl_id)]
    # chembl_id_rxn = chembl_id_rxn.reset_index(drop=True)
    # for i in range(len(chembl_id_rxn)):
        # itmpGraph.add_association(Abundance(namespace='ChEMBL',name = chembl_id_rxn['chembl_id'][i]),
                                  # Pathology(namespace='SideEffect',name = chembl_id_rxn['event'][i]),
                                  # citation = "OpenTargets Platform",evidence = 'DrugReactions')
        
    # return(itmpGraph)
                        
    
# def chembl2rxn_rel(itmpGraph):
    
    # infile = open('data/normalized_data/drugReactions.pkl','rb')
    # rxn_df = pickle.load(infile)
    # infile.close()
    
    # chembl_id = []
    # for node in itmpGraph.nodes():
        # if isinstance(node,pybel.dsl.Abundance):
            # if node.namespace == 'ChEMBL':
                # chembl_id.append(node.name)
            
    # chembl_id_rxn = rxn_df[rxn_df['chembl_id'].isin(chembl_id)]
    # chembl_id_rxn = chembl_id_rxn.reset_index(drop=True)
    # for i in range(len(chembl_id_rxn)):
        # itmpGraph.add_association(Abundance(namespace='ChEMBL',name = chembl_id_rxn['chembl_id'][i]),
                                  # Pathology(namespace='SideEffect',name = chembl_id_rxn['event'][i]),
                                  # citation = "OpenTargets Platform",evidence = 'DrugReactions')
        
    # return(itmpGraph)
                              

#import time
# def _get_target_data(
#     protein_list: list,
#     organism: str
# ) -> pd.DataFrame:
#     """Get chemical for target data from ChEMBL.

#     :param protein_list:
#     :param organism:
#     :return:
#     """
#     df_data = []

#     target = new_client.target
#     activity = new_client.activity

#     for protein in tqdm(protein_list, desc='Retrieving chemicals for proteins'):
#         if pd.isna(protein):
#             continue
#         try:
#             prot_data = [target.search(protein)[0]]

#             # Search for protein with same synonym
#             if prot_data == [None]:
#                 prot_data = target.filter(
#                     target_synonym__icontains=protein, target_organism__istartswith=organism
#                 ).only(['target_chembl_id', 'target_pref_name', 'molecule_chembl_id', 'molecule_pref_name'])
#         except HttpBadRequest:
#             print(f'Cannot search for {protein} due to chembl error')
#             continue

#         # No results found
#         if not prot_data:
#             continue

#         for prot in prot_data:
#             # Absence of chembl id
#             if not prot['target_chembl_id']:
#                 continue

#             prot_activity_data = activity.filter(
#                 target_chembl_id=prot['target_chembl_id'],
#                 assay_type_iregex='(B|F)',
#             ).only([
#                 'pchembl_value', 'molecule_chembl_id', 'activity_id', 'target_pref_name', 'molecule_pref_name'
#             ])

#             if len(prot_activity_data) < 1:
#                 continue

#             for i in prot_activity_data:
#                 tmp = {}

#                 try:
#                     if i['pchembl_value'] is None:
#                         continue
#                 except HttpApplicationError:
#                     continue

#                 pchembl_val = i['pchembl_value']

#                 if float(pchembl_val) < 6:
#                     tmp['activity'] = 'inhibitor'
#                 else:
#                     tmp['activity'] = 'activator'

#                 tmp.update({
#                     'protein_symbol': protein,
#                     'protein_name': i['target_pref_name'],
#                     'aid': str(i['activity_id']),
#                     'chembl_id': i['molecule_chembl_id'],
#                     'compound_name': i['molecule_pref_name'].capitalize() if i['molecule_pref_name'] else ''
#                 })
#                 df_data.append(tmp)
            
#             time.sleep(5)
#             print('5 seconds of break')

#     # Merge duplicated protein-chemical entries into one
#     df = pd.DataFrame()
#     if len(df_data)==0:
#         return df
        

#     for idx, row in tqdm(enumerate(df_data), total=len(df_data), desc='Preparing data'):
#         if idx == 0:
#             df = df.append(row, ignore_index=True)
#         else:
#             _in_df = df.loc[
#                 (df['protein_symbol'] == row['protein_symbol']) & (df['chembl_id'] == row['chembl_id'])
#                 ]

#             if _in_df.empty:
#                 df = df.append(row, ignore_index=True)
#             else:
#                 row_index = _in_df.index

#                 # Check existing citations
#                 existing_assays = set(df.loc[row_index, 'aid'].values[0].split(' | '))
#                 old_count = len(existing_assays)
#                 existing_assays.add(row['aid'])
#                 new_count = len(existing_assays)

#                 # Check if new citation added, if yes - add respective data
#                 if old_count < new_count:
#                     df.loc[row_index, 'aid'] = ' | '.join(existing_assays)
#     print(df)
#     print(df_data)
#     df = df[['activity', 'protein_symbol', 'protein_name', 'aid', 'chembl_id', 'compound_name']]
    
#     return df


# def target_list_to_chemical(
#     proteins: list,
#     organism: str = 'Homo sapiens',
# ) -> pd.DataFrame:
#     """Extract chemical information on list of targets
#     Usage:
#     >> target_list_to_chemical(proteins=['RIPK'])
#     """

#     df = _get_target_data(protein_list=proteins, organism=organism)
#     return df

#def _get_target_data(protein_list: list, organism: str):
#     """Get chemical for target data from ChEMBL"""
#     df_data = []

#     target = new_client.target
#     activity = new_client.activity

#     for protein in protein_list:
#         if pd.isna(protein):
#             continue
#         try:
#             prot_data = [target.search(protein)[0]]

#             # Search for protein with same synonym
#             if prot_data == [None]:
#                 prot_data = target.filter(
#                     target_synonym__icontains=protein, target_organism__istartswith=organism
#                 ).only(['target_chembl_id', 'target_pref_name', 'molecule_chembl_id', 'molecule_pref_name'])
#         except HttpBadRequest:
#             print(f'Cannot search for {protein} due to chembl error')
#             continue

#         # No results found
#         if not prot_data:
#             continue

#         for prot in tqdm(prot_data, f'Analying data for {protein}'):
#             # Absence of chembl id
#             if not prot['target_chembl_id']:
#                 continue

#             prot_activity_data = activity.filter(
#                 target_chembl_id=prot['target_chembl_id'],
#                 assay_type_iregex='(B|F)',
#             ).only([
#                 'pchembl_value', 'molecule_chembl_id', 'activity_id', 'target_pref_name', 'molecule_pref_name'
#             ])

#             if len(prot_activity_data) < 1:
#                 continue

#             for i in prot_activity_data:
#                 tmp = {}

#                 if i['pchembl_value'] is None:
#                     continue

#                 pchembl_val = i['pchembl_value']

#                 if float(pchembl_val) < 6:
#                     tmp['activity'] = 'inhibitor'
#                 else:
#                     tmp['activity'] = 'activator'

#                 tmp['protein_symbol'] = protein
#                 tmp['protein_name'] = i['target_pref_name']
#                 tmp['aid'] = str(i['activity_id'])
#                 tmp['chembl_id'] = i['molecule_chembl_id']
#                 tmp['compound_name'] = i['molecule_pref_name'].capitalize() if i['molecule_pref_name'] else ''
#                 df_data.append(tmp)
                
#             time.sleep(10)

#     # Merge duplicated protein-chemical entries into one
#     df = pd.DataFrame()

#     for idx, row in tqdm(enumerate(df_data), total=len(df_data), desc='Preparing data'):
#         if idx == 0:
#             df = df.append(row, ignore_index=True)
#         else:
#             _in_df = df.loc[
#                 (df['protein_symbol'] == row['protein_symbol']) & (df['chembl_id'] == row['chembl_id'])
#             ]

#             if _in_df.empty:
#                 df = df.append(row, ignore_index=True)
#             else:
#                 row_index = _in_df.index

#                 # Check existing citations
#                 existing_assays = set(df.loc[row_index, 'aid'].values[0].split(' | '))
#                 old_count = len(existing_assays)
#                 existing_assays.add(row['aid'])
#                 new_count = len(existing_assays)

#                 # Check if new citation added, if yes - add respective data
#                 if old_count < new_count:
#                     df.loc[row_index, 'aid'] = ' | '.join(existing_assays)
#     df = df[['activity', 'protein_symbol', 'protein_name', 'aid', 'chembl_id', 'compound_name']]
#     return df

# def target_list_to_chemical(
#     proteins: list,
#     organism: str = 'Homo sapiens',
#     output_dir: str = ''
# ) -> None:
#     """Extract chemical information on list of targets
#     Usage:
#     >> target_list_to_chemical(proteins=['RIPK'])
#     """

#     df = _get_target_data(protein_list=proteins, organism=organism)
#     #os.makedirs(output_dir, exist_ok=True)
#     #df.to_csv(os.path.join(output_dir, 'chemical_annotated.csv'), sep='\t', index=False)
#     return(df)
   