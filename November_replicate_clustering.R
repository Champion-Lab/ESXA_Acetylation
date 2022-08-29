library(flexclust)
library(factoextra)


#clustering with ratios
Nov_cluster_acet <- data.frame(protein = final_Nov$accession,
                               WT = final_Nov$WT_av_acet,
                               Del = final_Nov$Del_av_acet,
                               Comp = final_Nov$Comp_av_acet)
row.names(Nov_cluster_acet) <- Nov_cluster_acet$protein
Nov_cluster_acet <- subset(Nov_cluster_acet, select = -(protein))
Nov_cluster_acet <- Nov_cluster_acet[rowSums(Nov_cluster_acet) > 0,]

Nov_cluster_acet_temp <- Nov_cluster_acet

Nov_cluster_acet$Del_Comp_rat <- Nov_cluster_acet$Del / Nov_cluster_acet$Comp
Nov_cluster_acet$Del_WT_rat <- Nov_cluster_acet$Del / Nov_cluster_acet$WT
Nov_cluster_acet$Comp_WT_rat <- Nov_cluster_acet$Comp / Nov_cluster_acet$WT

Nov_cluster_acet <- subset(Nov_cluster_acet, select = -c(WT, Del, Comp))
Nov_cluster_acet <- as.data.frame(scale(Nov_cluster_acet, center = F))

#creates the clustered object
#Nov_cluster_acet is the name of the dataframe with data, 3 is number of clusters
# 25 is the number of starting positions (higher nstart takes longer but produces more
# reliable results)
Nov_kmeans_acet <- kmeans(Nov_cluster_acet, 3, nstart = 25)
#print the clustered object
Nov_kmeans_acet


#add the clusters to the original dataframe, save as new object
Nov_cluster_acet_4 <- cbind(Nov_cluster_acet_temp, cluster = Nov_kmeans_acet$cluster)

#visualize the clusters (basically creates a PCA plot)
clusterplot <- fviz_cluster(Nov_kmeans_acet, Nov_cluster_acet, ggtheme = theme_bw(), geom = "point",
                            ellipse = T, main = F) +
  theme_bw(base_size = 20) +
  scale_color_manual(values = c("darkgoldenrod1", "magenta", "blue")) +
  scale_fill_manual(values = c("darkgoldenrod1", "magenta", "blue"))
#view the plot
clusterplot

#everything before the theme_bw line
clusterplot_to_save <- clusterplot <- fviz_cluster(Nov_kmeans_acet, Nov_cluster_acet, ggtheme = theme_bw(), geom = "point",
                                                   ellipse = T, main = F)
#save(clusterplot_to_save, file = "exports_2022_05_23/clusterplot_ggplot_obj.rdata")
#load("exports_2022_05_23/clusterplot_ggplot_obj.rdata")


Nov_cluster_acet_4$protein <- row.names(Nov_cluster_acet_4)
Nov_cluster_acet_4_long <- pivot_longer(Nov_cluster_acet_4, c(WT, Del, Comp), names_to = "condition", values_to = "area")
Nov_cluster_acet_4_long$cluster <- as.factor(Nov_cluster_acet_4_long$cluster)

Nov_cluster_acet_4_long$area[Nov_cluster_acet_4_long$area == 0] <- 0.00001

plot <- ggplot(Nov_cluster_acet_4_long) +
  geom_point(aes(x = condition, y = area, color = cluster, group = protein))+
  geom_line(aes(x = condition, y = area, color = cluster, group = protein), lwd = 1) +
  scale_y_continuous(trans = "log2",
                     labels = function(x) format(x, scientific = T),
                     breaks = c(1e-05, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 4)) +
  scale_color_manual(values = c("darkgoldenrod1", "magenta", "blue")) +
  theme_bw(base_size = 20) +
  labs(y = "Normalized Area", x = element_blank()) +
  theme(panel.grid = element_blank(), legend.position = "none")
plot


write.csv(Nov_cluster_acet_4_long, "exports_2022_05_23/plotexports_2022_05_31.csv", row.names = F)
cluster1df <- Nov_cluster_acet_4_long[Nov_cluster_acet_4_long$cluster == 1,]
cluster2df <- Nov_cluster_acet_4_long[Nov_cluster_acet_4_long$cluster == 2,]
cluster3df <- Nov_cluster_acet_4_long[Nov_cluster_acet_4_long$cluster == 3,]

plot1 <- ggplot(cluster1df) +
  geom_point(aes(x = condition, y = area, group = protein), color = "darkgoldenrod1")+
  geom_line(aes(x = condition, y = area, group = protein), lwd = 1, color = "darkgoldenrod1") +
  scale_y_continuous(trans = "log2",
                     labels = function(x) format(x, scientific = T),
                     breaks = c(1e-05, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 4),
                     limits = c(1e-05, 4)) +
  theme_bw(base_size = 20) +
  labs(y = "Normalized Area", x = element_blank()) +
  theme(panel.grid = element_blank())
plot1

plot2 <- ggplot(cluster2df) +
  geom_point(aes(x = condition, y = area, group = protein), color = "magenta")+
  geom_line(aes(x = condition, y = area, group = protein), lwd = 1, color = "magenta") +
  scale_y_continuous(trans = "log2",
                     labels = function(x) format(x, scientific = T),
                     breaks = c(1e-05, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 4),
                     limits = c(1e-05, 4)) +
  theme_bw(base_size = 20) +
  labs(y = "Normalized Area", x = element_blank()) +
  theme(panel.grid = element_blank())
plot2

plot3 <- ggplot(cluster3df) +
  geom_point(aes(x = condition, y = area, group = protein), color = "blue")+
  geom_line(aes(x = condition, y = area, group = protein), lwd = 1, color = "blue") +
  scale_y_continuous(trans = "log2",
                     labels = function(x) format(x, scientific = T),
                     breaks = c(1e-05, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 4),
                     limits = c(1e-05, 4)) +
  theme_bw(base_size = 20) +
  labs(y = "Normalized Area", x = element_blank()) +
  theme(panel.grid = element_blank())
plot3

