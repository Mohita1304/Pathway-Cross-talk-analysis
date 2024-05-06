fh1=open("Selected_Pathway_KEGG_geneset.csv","r")
pathway=fh1.readlines()
# print(pathway)

fh3=open("Type1_PPC_Tumor.txt","w")
# type1_crosstalk=[]
fh2=open("Gene_Coexp_Net.csv","r")
for line in fh2:
    line1=line.rstrip().split(" ")
    # print(line1)
    # print(line1[2])
    pathA=[]
    pathB=[]
    for path in pathway:
        Path_geneset=path.rstrip().split(",")
        # print(Path_geneset)
        if (line1[0] in Path_geneset) and (line1[1] not in Path_geneset):
            # print(line1[0])
            pathA.append(Path_geneset[0])
            # print("GeneI:",line1[0], Path_geneset) 
        if (line1[1] in Path_geneset) and (line1[0] not in Path_geneset):
            pathB.append(Path_geneset[0])
            # print(pathB)
    if ((len(pathA)!= 0) and (len(pathB)!=0)):
            for i in pathA:
                #print(i)
                for j in  pathB:
                    if i != j:
                            print(line1[0],i,line1[1], j)
                            fh3.write(line1[0]+ "::" + line1[1]+"\t"+ str(i)+ "\t"+str(j)+"\n")
            
       
        
        



