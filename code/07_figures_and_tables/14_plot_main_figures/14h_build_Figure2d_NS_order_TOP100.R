options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))
suppressPackageStartupMessages({library(ggplot2);library(scales)})
ROOT <- "__WGPE_PROJECT_ROOT__"
OUT <- file.path(ROOT, "output_TOP100_revision", "figures", "main",
                 "Figure2d_TOP100_NS_ORDER")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(file.path(ROOT, "output_TOP100_revision", "derived_data",
                        "figure_source_data",
                        "Figure2_ENSO_standard_band_decomposition_TOP100.csv"),
              check.names = FALSE)
display <- c("Global", "N high latitude", "N temperate", "Tropics",
             "S temperate", "S high latitude")
z <- rbind(
  data.frame(region=d$region,component="IWV",value=d$IWV_channel_J_m2),
  data.frame(region=d$region,component="<z>",value=d$z_channel_J_m2),
  data.frame(region=d$region,component="Interaction",value=d$interaction_channel_J_m2)
)
z$region <- factor(z$region, levels=rev(display))
z$component <- factor(z$component, levels=c("IWV","<z>","Interaction"))
p <- ggplot(z,aes(value,region,fill=component))+
  geom_vline(xintercept=0,colour="#606060",linewidth=.28)+
  geom_col(width=.52,colour="#222222",linewidth=.06)+
  scale_fill_manual(values=c("IWV"="#2C7BB6","<z>"="#D66A55","Interaction"="#8E2A62"),name=NULL)+
  scale_x_continuous(breaks=c(-5000,0,10000,20000,30000),labels=label_number(big.mark=" "))+
  coord_cartesian(xlim=c(-8000,36000),clip="on")+
  labs(title=expression(bold("d")~"ENSO decomposition"),
       x=expression("Channel anomaly (J m"^{-2}*")"),y=NULL)+
  theme_classic(base_size=9,base_family="sans")+
  theme(panel.grid.major.x=element_line(colour="#ECECEA",linewidth=.18),
        panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),
        axis.line.y=element_blank(),axis.ticks.y=element_blank(),
        legend.position="bottom",legend.direction="horizontal",
        plot.title=element_text(size=10.5,margin=margin(b=2)),
        plot.margin=margin(4,7,4,5))
ggsave(file.path(OUT,"Figure2d_ENSO_decomposition_NS_order_TOP100.png"),p,
       width=92,height=66,units="mm",dpi=600,bg="white")
writeLines(c("# Figure 2d north-to-south ordering QA","",
             "- Scientific source changed: NO.",
             "- Display order: Global, N high latitude, N temperate, Tropics, S temperate, S high latitude.",
             "- ONI-standard TOP100 exact channel decomposition retained: YES.",
             "- Output: PNG only."),
           file.path(OUT,"Figure2d_NS_order_QA.md"),useBytes=TRUE)
