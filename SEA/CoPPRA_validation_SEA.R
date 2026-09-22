# Plot helpers for sections 7.1 and 7.2. Typography matches Base QC (14/16/18 pt).
validation_theme_SEA <- function() {
  ggplot2::theme_set(ggplot2::theme_bw(base_size = 14) + ggplot2::theme(
    plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_blank(),
    plot.caption = ggplot2::element_blank(),
    plot.tag = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = 16), axis.text = ggplot2::element_text(size = 14),
    strip.text = ggplot2::element_text(size = 16, face = "bold")))
}
split_ids_validation <- function(x) {
  z <- trimws(unlist(strsplit(as.character(x[!is.na(x)]), ";", fixed = TRUE), use.names = FALSE))
  sort(unique(z[nzchar(z)]))
}
mod_categories <- function(x) {
  if (is.na(x) || !nzchar(x)) return("Unknown")
  sort(unique(trimws(sub("^[0-9]+ +", "", strsplit(x, ",", fixed = TRUE)[[1]]))))
}
canonical_group <- function(x) vapply(x, function(z) paste(split_ids_validation(z), collapse = ";"), character(1))
venn_validation <- function(sets, title) {
  stopifnot(identical(names(sets), c("Base", "SEA")))
  angle <- seq(0, 2*pi, length.out = 361)
  circles <- rbindlist(lapply(seq_along(sets), function(i)
    data.table(x = c(-0.65, 0.65)[i] + 1.2*cos(angle), y = 1.2*sin(angle), Search = names(sets)[i])))
  counts <- data.table(x = c(-1.1,0,1.1), y = 0,
    N = c(length(setdiff(sets$Base, sets$SEA)),length(intersect(sets$Base,sets$SEA)),length(setdiff(sets$SEA,sets$Base))))
  labs <- data.table(x = c(-1,1), y=1.5, Label=display_label_SEA(names(sets)))
  ggplot2::ggplot(circles, ggplot2::aes(x,y,group=Search,fill=Search,colour=Search)) +
    ggplot2::geom_polygon(alpha=.25,linewidth=.9) +
    ggplot2::geom_text(data=counts,ggplot2::aes(x,y,label=format(N,big.mark=",")),inherit.aes=FALSE,size=6) +
    ggplot2::geom_text(data=labs,ggplot2::aes(x,y,label=Label),inherit.aes=FALSE,size=5) +
    ggplot2::scale_fill_manual(values=c(Base="#2878B5",SEA="#E88931")) +
    ggplot2::scale_colour_manual(values=c(Base="#2878B5",SEA="#E88931")) +
    ggplot2::coord_fixed(xlim=c(-2.1,2.1),ylim=c(-1.35,1.85),clip="off") +
    ggplot2::labs(title=display_label_SEA(title),x=NULL,y=NULL) +
    ggplot2::theme(legend.position="none",axis.text=ggplot2::element_blank(),
      axis.ticks=ggplot2::element_blank(),panel.grid=ggplot2::element_blank(),
      plot.subtitle=ggplot2::element_blank(),plot.caption=ggplot2::element_blank(),
      plot.tag=ggplot2::element_blank())
}
mod_summary_validation <- function(ev, dataset, scope) {
  forms <- unique(ev[, .(Sequence, Modified.sequence, Modifications)])
  label_map <- rbindlist(lapply(unique(forms$Modifications), function(z)
    data.table(Modifications=z,Category=mod_categories(z))))
  expanded <- merge(forms,label_map,by="Modifications",allow.cartesian=TRUE)
  out <- expanded[, .(Sequences=uniqueN(Sequence),Modified_forms=uniqueN(Modified.sequence)),by=Category]
  rows <- merge(ev[, .(Evidence_rows=.N),by=Modifications],label_map,by="Modifications",allow.cartesian=TRUE)[
    ,.(Evidence_rows=sum(Evidence_rows)),by=Category]
  out <- merge(out,rows,by="Category")
  out[, `:=`(Dataset=dataset,Scope=scope)]
  out[]
}
