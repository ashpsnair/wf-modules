

# Assuming your Seurat object is named 'seurat_obj'

mean_expression <- mean(GetAssayData(seurat_obj, slot = "data")["AMPK",])
upregulated_cells <- WhichCells(object = seurat_obj, expression = AMPK > mean_expression)

# Subset the Seurat object to only include cells with upregulated AMPK
seurat_obj_upregulated <- subset(seurat_obj, cells = upregulated_cells)

'''
split the bam based on the expression level
then further split based on the cell  type to run the scomatic pipelinwe

'''