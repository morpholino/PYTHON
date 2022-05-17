#make a pallette of colours for a set of locations
from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import csv
import os
import argparse

#################
### FUNCTIONS ###
#################


def my_branch_style():
    othersupport, fullsupport, highsupport, leaf = NodeStyle(), NodeStyle(), NodeStyle(), NodeStyle()
    fullsupport["hz_line_width"] = 5
    fullsupport["size"] = 0
    #fullsupport["extra_branch_line_type"] = 2 #the extra branches need to be adjusted still

    highsupport["hz_line_width"] = 5
    highsupport["hz_line_color"] = "#666666"
    highsupport["size"] = 0
    #highsupport["extra_branch_line_type"] = 2 #the extra branches need to be adjusted still
    othersupport["size"] = 0

    return  othersupport, fullsupport, highsupport, leaf


def my_tree_style(filename, circular):
    ts = TreeStyle()
    
    if circular:
        ts.mode = "c" # draw tree in circular mode
        ts.scale = 100
        ts.arc_start = -90 # 0 degrees = 3 o'clock
        ts.arc_span = 330
    else:
        ts.show_leaf_name = True
        ts.min_leaf_separation = 1
        #ts.show_branch_length = True
        ts.show_branch_support = True
        ts.scale =  200 # 100 pixels per branch length unit
        ts.min_leaf_separation = 0.5
        ts.branch_vertical_margin = 0 # 10 pixels between adjacent branches
        ts.title.add_face(TextFace(filename, fsize=20), column=0)

    return ts


def find_supported(tree, support):
    #"Find nodes with significant support"
    matches = []
    for n in tree.traverse():
        if n.support >= support and n.is_leaf() == False and n.is_root() == False: 
            matches.append(n)
            n.add_features(name="SUP")
    return matches


def tag_replace(string):
    #print("Lineage lookup began, please use the following list of unrecognized genera to update your fetch_lineages.tsv")

    if "-" in string:
        string = string.replace("-", "_")
    if " " in string:
        string = string.replace(" ", "_")
    tag = string.split("_")[0]
    if tag in {"Candidatus", "uncultured", "mt", "pt", "cyto", "SEED", "SEED1", "SEED2", "WEIRD", "candidate"}:
        tag = string.split("_")[1]
    if tag in taxarepl9:
        genus = taxarepl9.get(tag, "").split(" ")[0]
        string = string.replace(tag, taxarepl9[tag])
    else:
        genus = tag
    if tag in fetch_lineages:
        string = "{}@{}".format(string, fetch_lineages[tag])
    elif genus in fetch_lineages:
        string = "{}@{}".format(string, fetch_lineages[genus])
    else:
        print("No lineage available for: {}".format(tag))
        all_lineages_found = False

    return string
   

#################
### DATA READ ###
#################


def main():
    parser = argparse.ArgumentParser(description='How to use argparse')
    parser.add_argument('-t', '--treefile', help='Tree file to be processed', default="")
    parser.add_argument('-a', '--ancestors', help='Ancestors file in tsv format', default="_ancestors.tsv")
    parser.add_argument('-l', '--localization', help='Localization data', default="")
    parser.add_argument('-c', '--circular', help='Circular tree output', action="store_true")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-s', '--some_other_function', help='Maybe something else with the tree', default=None)
    group.add_argument('-p', '--prune', help='Prune mask', default=True)

    args = parser.parse_args()

    if args.treefile == "":
        filename = "_matrix_metamonada"
    else:
        filename = args.treefile

    if args.prune:
        prunefile = "_matrix_metamonada_skip_cafe.txt"
        #prunefile = args.prune

    if os.path.isdir("/Users/morpholino/OwnCloud/"):
        home = "/Users/morpholino/OwnCloud/"
    elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
        home = "/Volumes/zoliq data/OwnCloud/"
    else:
        print("Please set a homedir")
    if True:
        print("setting default directory")
        coredata = home + "progs/PYTHON/"
        wd = home + "progs/PYTHON/ETEtree/"

    #load additional data
    taxareplfile = csv.reader(open(coredata + "taxarepl9.tsv"), delimiter="\t", skipinitialspace=False)
    taxarepl9 = {row[0]:row[1] for row in taxareplfile}
    lineagefile = csv.reader(open(coredata + "fetch_lineages.tsv"), delimiter="\t", skipinitialspace=False)
    fetch_lineages = {row[0]:row[1] for row in lineagefile}

    all_lineages_found = True
    #################
    print("loading tree file")
    t = Tree(filename + ".treefile", format=0) #format flexible with support values


    ##################
    ###    MAIN    ###
    ##################
    #set the ancestors here...
    ancestors_d = {"Clade1": ["cyanANABv_WP_011317903.1_plastoquinolTOX", 
                              "grnPRASc_MMETSP0941_Gene.14464-Transcript_5625_Chlorophyta_Prasinococcales"],
                   "Clade2": ["crypGONIp_MMETSP0107_Gene.30083-Transcript_20766_Cryptophyta_Cryptomonadales",
                              "Phenylobacterium_sp._RIFCSPHIGHO2_01_FULL_69_31_OHB27812.1"]
                  }
    #or make them available in a tsv:
    with open(args.ancestors, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        try:
            ancestors_d.update({r[0]: [r[1], r[2]] for r in reader})
        except IndexError:
            #incomplete line
            print("Incomplete processing of ancestor lineages!")
            print(ancestors_d)

    try:
        ancestor = t.get_common_ancestor(ancestors_d[filename])
        #print(ancestor)
        t.set_outgroup(ancestor)
    except KeyError:
        print("Root not selected!")
        print(t.get_tree_root())
        quit()
    t.ladderize(direction=1)


    #################################
    ###    LOCALIZE SUBROUTINE    ###
    #################################
     #create a dictionary of taxon data and a list of all localities
    if args.localization:
        print("loading data files")  
        LocDataFile = csv.reader(open(wd + "4pred-preds.txt"), delimiter="\t", skipinitialspace=True)

        locdata = {}
        localization = set()

        for row in LocDataFile:
            locdata[row[0]] = row[2]
            localization.add(row[2])

        list(localization).sort()

        #assign colors to localizations
        N = len(localization)
        colors = {"MT": "dodgerblue", "PT": "mediumseagreen", "dual": "darkviolet", 
                  "CS": "black", "amb": "black", 
                  "SP": "black", "no pred": "black"}
        #automatic assignment:
        #cmap = get_cmap(N)
        #colors = {}
        #for i, X in enumerate(localization):
        #    colors[X] = cmap[i] #hopefully this is still recognized as a color
        #print(colors)


    ################################
    ###    SAVEFIG SUBROUTINE    ###
    ################################
    #select scale 0-1.0 or 0-100 for support values
    supportscache = t.get_cached_content(store_attr="support")
    supportslist = [x.support for x in supportscache]
    if max(supportslist) == 1:
        minsupport = 0.85
    else:
        minsupport = 85
    find_supported(t, support=minsupport) #find non-terminal nodes with high support


    #set different graphic styles:
    othersupport, fullsupport, highsupport, leaf = my_branch_style()

    #define tree style
    ts = my_tree_style(filename, args.circular)

    pruned = []
    #now for something completely different, tree traverse
    for node in t.traverse():
        if node.is_leaf():
            pruned.append(node)
            if args.localization:
                local = locdata.get(node.name, "no pred")
                leafcolor = colors.get(local, "black")
                node.add_features(color=leafcolor)
                node.img_style["fgcolor"] = leafcolor
                node.img_style["vt_line_color"] = leafcolor
                node.img_style["hz_line_color"] = leafcolor
            node.img_style["size"] = 0
            node.img_style["hz_line_width"] = 2
            #node.name = tag_replace(node.name) #make s
            #print(node.name)
            #node.set_style(leaf)
            #print(node.name, local)
        elif node.name == "SUP":
            if node.support == 100:
                node.set_style(fullsupport)
            else:
                node.set_style(highsupport)
        else:
            node.set_style(othersupport)
    #to check tree is annotated:
    #print(t.get_ascii(attributes=["name", "abundance"], show_internal=True)) 
    #color and abundance are new features

    print("all terminal branches:", len(pruned))
    if args.prune:
        with open(prunefile) as f:
            skipnodes = f.read().split("\n")
            skipnodes = [x.split("[&!")[0] for x in skipnodes]
            print("To be pruned:\n{}\n\n".format(len(skipnodes), "\n".join(skipnodes)))

        pruned = [x for x in pruned if x.name not in skipnodes]
        print("remain after pruning:", len(pruned))
        t.prune(pruned, preserve_branch_length=True)

    #ts.legend.add_face(fig, column=0) #does not work
    #t.show(tree_style=ts)
    t.write(format=0, outfile="{}_pruned.treefile".format(filename))
    t.render(filename + "_c.pdf", w=400, units="mm", tree_style=ts)

    if all_lineages_found == False:
        print("run fetch_entrez_lineage_tree.py to find missing lineages")

    print("Analysis finished")


if __name__ == "__main__":
    main()