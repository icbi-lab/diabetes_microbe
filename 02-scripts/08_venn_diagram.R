library(VennDiagram)


T1DM_vs_H <- c("Blautia", "Subdoligranulum", "CAG-352", "Akkermansia", "Fusicatenibacter")
T3cDM_vs_H <- c("Subdoligranulum", "Streptococcus", "Akkermansia", "Blautia", "CAG-352")
T3cDM_vs_T1DM <- c("Streptococcus", "Blautia", "Escherichia-Shigella", "Faecalibacterium", "CAG-352")

venn_list <- list(
  "T1DM vs H" = T1DM_vs_H,
  "T3cDM vs H" = T3cDM_vs_H,
  "T3cDM vs T1DM" = T3cDM_vs_T1DM
)

venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = c("#A3D1D1", "#FFFF66", "#FFD1D1"),
  alpha = 0.6,
  col = "black",
  cex = 2.2,          
  cat.cex = 1.8,          
  cat.col = "black",
  cat.pos = c(-20, 20, 0),
  cat.dist = c(0.05, 0.05, 0.05),
  margin = 0.1,
  main = ""
)

# Draw to screen
grid.newpage()
grid.draw(venn.plot)

#ggsave(plot=p,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/venn_diagram.svg", height = 8, width = 15)
#ggsave(plot=p,"/data/scratch/kvalem/projects/2024/diabetes_microbe/05-results/figures/venn_diagram.png", height = 8, width = 15)

