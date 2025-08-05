library(tidyverse)
library(dplyr)

length <- read_tsv("all.contigs.Length.txt") %>% group_by(length) %>% summarise(Count = n())
sum <- sum(length$Count)

interval <- length
# 创建一个新的数据框，其中包含分组变量和count的总和
# grouped_df <- interval %>%
#   mutate(
#     # 使用case_when创建分组变量
#     length_group = case_when(
#       length < 200 ~ "0-199",
#       length >= 200 & length < 2000 ~ "200-1999",
#       length >= 2000 & length < 5000 ~ "2000-5000",
#       length >= 5000 ~ "5000+",
#     )
#   ) %>%
#   group_by(length_group) %>%
#   summarise(
#     total_count = sum(Count),
#     .groups = 'drop' # 删除分组属性，以便结果是一个普通的数据框
#   )

grouped_df <- interval %>%
  mutate(
    # 使用case_when创建分组变量
    length_group = case_when(
      length < 200 ~ "0-199",
      length >= 200 & length < 300 ~ "200-299",
      length >= 300 & length < 400 ~ "300-399",
      length >= 400 & length < 500 ~ "400-499",
      length >= 500 & length < 600 ~ "500-599",
      length >= 600 & length < 700 ~ "600-699",
      length >= 700 & length < 800 ~ "700-799",
      length >= 800 & length < 900 ~ "800-899",
      length >= 900 & length < 1000 ~ "900-999",
      length >=1000 & length < 1100 ~ "1000-1099",
      length >= 1100 & length < 1200 ~ "1100-1199",
      length >= 1200 & length < 1300 ~ "1200-1299",
      length >= 1300 & length < 1400 ~ "1300-1399",
      length >= 1400 & length < 1500 ~ "1400-1499",
      length >= 1500 & length < 1600 ~ "1500-1599",
      length >= 1600 & length < 1700 ~ "1600-1699",
      length >= 1700 & length < 1800 ~ "1700-1799",
      length >= 1800 & length < 1900 ~ "1800-1899",
      length >= 1900 & length <= 2000 ~ "1900-1999",
      length >= 2000 & length < 5000 ~ "2000-5000",
      length >= 5000 ~ "5000+",
    )
  ) %>%
  group_by(length_group) %>%
  summarise(
    total_count = sum(Count),
    .groups = 'drop'# 删除分组属性，以便结果是一个普通的数据框
  )


# 计算整体的总和
total_overall <- sum(interval$Count)

# 添加占比列
grouped_df <- grouped_df %>%
  mutate(
    percentage = (total_count / total_overall) * 100
  )

grouped_df$percentage <- round(grouped_df$percentage)


p <- ggplot(length) + 
  geom_col(aes(length, Count), width = 0.8) + 
#  geom_line(aes(length, Count), group = 1) + 
#  geom_point(aes(length, Count)) + 
  scale_y_continuous(sec.axis = sec_axis(~.*100/sum, name = "% Relative Abundance"),expand=c(0, 0)) + 
  xlab("Length") +
  theme_bw() + theme(panel.grid = element_blank(), 
                     axis.title = element_text(size = 15))+
  scale_x_continuous(limits = c(200,2000), breaks = seq(200, 2000, 100),expand=c(0, 0))

p
ggsave(p,filename = "osahs.all.contigs.Length200-2000.pdf", height = 5, width = 8)



grouped_df$length_group<- factor(grouped_df$length_group,levels =c('0-199',
                                                                   '200-299',
                                                                   '300-399',
                                                                   '400-499',
                                                                   '500-599',
                                                                   '600-699',
                                                                   '700-799',
                                                                   '800-899',
                                                                   '900-999',
                                                                   '1000-1099',
                                                                   '1100-1199',
                                                                   '1200-1299',
                                                                   '1300-1399',
                                                                   '1400-1499',
                                                                   '1500-1599',
                                                                   '1600-1699',
                                                                   '1700-1799',
                                                                   '1800-1899',
                                                                   '1900-1999',
                                                                   '2000-5000',
                                                                   '5000+'))
p2 <- ggplot(grouped_df, aes(x = length_group, y = percentage)) +      
  geom_bar(stat = "identity",width = 0.6)+
  theme_classic()+
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust = 0.5, vjust = 0.5))+ geom_text(aes(label = paste(percentage,'%')), vjust = -0.5,
                                                         color = "black", size = 3.5) +
  scale_y_continuous(limits = c(0,100),expand=c(0, 0))

p2
ggsave(p2,filename = "osahs.all.contigs.Length_group_2.pdf", height = 5, width = 8)
