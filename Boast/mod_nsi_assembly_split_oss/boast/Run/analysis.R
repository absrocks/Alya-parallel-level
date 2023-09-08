#   Copyright (C) 2017  Damien Dosimont
#
#   This a is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(ggplot2)
library(RColorBrewer)


cheader<-c("type", "vector", "inline", "unroll", "optimization", "time")


read <- function(file, cheader, sep=',') {
  df <- read.csv(file, header=TRUE, sep = sep, strip.white=TRUE)
  names(df) <- cheader
  df
}

comparison <- function(data){
  data$label <- paste(data$vector, data$inline, data$unroll, data$optimization, sep=" ")
  plot <- ggplot(data, aes(label, time, fill=type)) + theme(axis.text.x = element_text(angle=60, hjust=1))
  plot <- plot + geom_bar(stat="identity", position="dodge")
  plot <- plot + ggtitle("Time vs Parameters")
  plot <- plot + labs(x="Parameters", y="Time (seconds)", fill="Boast/Ref")
  plot
}

select_best <- function(data){
  #data <- data[data$type %in% "boast",]
  data$label <- paste(data$type, data$vector, data$inline, data$unroll, data$optimization, sep=" ")
  data$label <- factor(data$label, levels = unique(data$label)[order(data$time)])
  plot <- ggplot(data, aes(label, time, fill=type)) + theme(axis.text.x = element_text(angle=60, hjust=1))
  plot <- plot + geom_bar(stat="identity", position="dodge")
  plot <- plot + ggtitle("Time vs Partition/IO")
  plot <- plot + labs(x="Parameters", y="Time (seconds)", fill="Boast/Ref")
  plot
}

#managing the arguments
args <- commandArgs(trailingOnly = TRUE)

#file
file <- "test.csv"

#data
data <-read(file, cheader, ',')

comparison(data)
ggsave("res_comparison.pdf")
select_best(data)
ggsave("res_sorted.pdf")

#warnings()
