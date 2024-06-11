# Tyep II pathway cross-talk

fh1=open("Selected_Pathway_KEGG_geneset.csv","r")
pathway=fh1.readlines()
#print
fh1.close()

fh2=open("Gene_Coexpression_input.csv","r") 
GenePair=fh2.readlines()
# print(GenePair)
fh2.close()


fh3=open("Type2_PPC_Tumor.csv", "w")
for k, line in enumerate(GenePair):
    Gene=line.rstrip().split(" ")
    # print(Gene)
    GeneA=set(Gene[0:2])
    # print(GeneA)
    for line1 in GenePair[k:]:
        Gene1=line1.rstrip().split(" ")
        # print(Gene1)
        GeneB=set(Gene1[0:2])
        # print(GeneB)
        if Gene[0] != Gene1[0] and Gene[1] != Gene1[1]:
            # print(Gene,",", Gene1)
            shared_Gene = GeneA & GeneB
            # print(shared_Gene)
            sharedGene = list(GeneA & GeneB)
            # print(sharedGene)
            Shared_Path = []
            Uniq_PathA = []
            Uniq_PathB = []
            if len(shared_Gene) > 0:
                # print(shared_Gene)
                uniq_GeneA = list(GeneA - GeneB)
                # print(uniq_GeneA)
                uniq_GeneB = list(GeneB - GeneA)
                # print(uniq_GeneB)
                for path in pathway:
                    Path_geneset = path.rstrip().split(",")
                    # print(Path_geneset)
                    if ((uniq_GeneA[0] in Path_geneset) and (sharedGene[0] in Path_geneset) and (uniq_GeneB[0] not in Path_geneset)):
                        Uniq_PathA.append(Path_geneset[0])
                        #print(Uniq_PathA)
                    if ((uniq_GeneB[0] in Path_geneset) and (sharedGene[0] in Path_geneset) and (uniq_GeneA[0] not in Path_geneset) ):
                        Uniq_PathB.append(Path_geneset[0])
                    # print(Uniq_PathB)
                for i in Uniq_PathA:
                    for j in Uniq_PathB:
                        if i != j:
                            print("GeneA:", GeneA, uniq_GeneA, i, "SharedGEne", sharedGene, "GeneB", GeneB, uniq_GeneB,j)
                            fh3.write('\t'.join(uniq_GeneA) + ";" + i + ";" + '\t'.join(sharedGene) + ";" + '\t'.join(
                                uniq_GeneB) + ";" + j + "\n")
fh3.close()
