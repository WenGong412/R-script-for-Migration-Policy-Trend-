# R-script-for-Migration-Policy-Trend-Wen-Gong
#rvest
install.packages('rvest')
library('rvest')

#devtools/quanteda
install.packages("devtools")
library(devtools)
library(quanteda)

#installing qdap&rJava
install_github("trinker/qdapDictionaries")
install_github("trinker/qdapRegex")
install_github("trinker/qdapTools")
install_github("trinker/qdap")
install.packages("rJava")
install.packages("qdap")
library(qdap)

url <- 'http://www.migrationpolicy.org/article/trump-administration-six-months-sea-change-immigration-enforcement'
webpage <- read_html(url)
description_data_html <- html_nodes(webpage,'p')
description_data <- html_text(description_data_html)
mydfm <- dfm(description_data, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

# url1
url1 <- 'http://www.migrationpolicy.org/article/philippines-beyond-labor-migration-toward-development-and-possibly-return'
webpage1 <- read_html(url1)
description_data_html1 <- html_nodes(webpage1,'p')
description_data1 <- html_text(description_data_html1)
mydfm1 <- dfm(description_data1, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url2 <- 'http://www.migrationpolicy.org/article/cuban-migration-postrevolution-exodus-ebbs-and-flows'
webpage2 <- read_html(url2)
description_data_html2 <- html_nodes(webpage2,'p')
description_data2 <- html_text(description_data_html2)
mydfm2 <- dfm(description_data2, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url3 <- 'http://www.migrationpolicy.org/article/immigrant-health-care-workers-united-states'
webpage3 <- read_html(url3)
description_data_html3 <- html_nodes(webpage3,'p')
description_data3 <- html_text(description_data_html3)
mydfm3 <- dfm(description_data3, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url4 <- 'http://www.migrationpolicy.org/article/despite-political-resistance-use-temporary-worker-visas-rises-us-labor-market-tightens'
webpage4 <- read_html(url4)
description_data_html4 <- html_nodes(webpage4,'p')
description_data4 <- html_text(description_data_html4)
mydfm4 <- dfm(description_data4, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url5 <- 'http://www.migrationpolicy.org/article/amid-economic-crisis-and-political-turmoil-venezuelans-form-new-exodus'
webpage5 <- read_html(url5)
description_data_html5 <- html_nodes(webpage5,'p')
description_data5 <- html_text(description_data_html5)
mydfm5 <- dfm(description_data5, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url6 <- 'http://www.migrationpolicy.org/article/refugees-and-asylees-united-states'
webpage6 <- read_html(url6)
description_data_html6 <- html_nodes(webpage6,'p')
description_data6 <- html_text(description_data_html6)
mydfm6 <- dfm(description_data6, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url7 <- 'http://www.migrationpolicy.org/article/estonian-citizenship-policy-restoration-country-leads-statelessness-some'
webpage7 <- read_html(url7)
description_data_html7 <- html_nodes(webpage7,'p')
description_data7 <- html_text(description_data_html7)
mydfm7 <- dfm(description_data7, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url8 <- 'http://www.migrationpolicy.org/article/texas-leads-resurgence-restrictive-state-actions-immigration-enforcement'
webpage8 <- read_html(url8)
description_data_html8 <- html_nodes(webpage8,'p')
description_data8 <- html_text(description_data_html8)
mydfm8 <- dfm(description_data8, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url9 <- 'http://www.migrationpolicy.org/article/russia-migration-system-soviet-roots'
webpage9 <- read_html(url9)
description_data_html9 <- html_nodes(webpage9,'p')
description_data9 <- html_text(description_data_html9)
mydfm9 <- dfm(description_data9, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url10 <- 'http://www.migrationpolicy.org/article/strife-abroad-responses-home-muslims-west-and-conflict-spillover'
webpage10 <- read_html(url10)
description_data_html10 <- html_nodes(webpage10,'p')
description_data10 <- html_text(description_data_html10)
mydfm10 <- dfm(description_data10, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url11 <- 'http://www.migrationpolicy.org/article/sub-saharan-african-immigrants-united-states'
webpage11 <- read_html(url11)
description_data_html11 <- html_nodes(webpage11,'p')
description_data11 <- html_text(description_data_html11)
mydfm11 <- dfm(description_data11, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url12 <- 'http://www.migrationpolicy.org/article/uptick-northern-border-crossings-places-canada-us-safe-third-country-agreement-under'
webpage12 <- read_html(url12)
description_data_html12 <- html_nodes(webpage12,'p')
description_data12 <- html_text(description_data_html12)
mydfm12 <- dfm(description_data12, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url13 <- 'http://www.migrationpolicy.org/article/uptick-northern-border-crossings-places-canada-us-safe-third-country-agreement-under'
webpage13 <- read_html(url13)
description_data_html13 <- html_nodes(webpage13,'p')
description_data13 <- html_text(description_data_html13)
mydfm13 <- dfm(description_data13, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url14 <- 'http://www.migrationpolicy.org/article/inmigrantes-centroamericanos-en-los-estados-unidos'
webpage14 <- read_html(url14)
description_data_html14 <- html_nodes(webpage14,'p')
description_data14 <- html_text(description_data_html14)
mydfm14 <- dfm(description_data14, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url15 <- 'http://www.migrationpolicy.org/article/despite-little-action-yet-trump-administration-sanctuary-cities-states-and-localities-rush'
webpage15 <- read_html(url15)
description_data_html15 <- html_nodes(webpage15,'p')
description_data15 <- html_text(description_data_html15)
mydfm15 <- dfm(description_data15, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url16 <- 'http://www.migrationpolicy.org/article/colombia-emerges-decades-war-migration-challenges-mount'
webpage16 <- read_html(url16)
description_data_html16 <- html_nodes(webpage16,'p')
description_data16 <- html_text(description_data_html16)
mydfm16 <- dfm(description_data16, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url17 <- 'http://www.migrationpolicy.org/article/central-american-immigrants-united-states'
webpage17 <- read_html(url17)
description_data_html17 <- html_nodes(webpage17,'p')
description_data17 <- html_text(description_data_html17)
mydfm17 <- dfm(description_data17, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url18 <- 'http://www.migrationpolicy.org/article/its-population-ages-japan-quietly-turns-immigration'
webpage18 <- read_html(url18)
description_data_html18 <- html_nodes(webpage18,'p')
description_data18 <- html_text(description_data_html18)
mydfm18 <- dfm(description_data18, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url19 <- 'http://www.migrationpolicy.org/article/muscular-public-relations-strategy-paint-immigrants-and-immigration-negatives-embedded-deep'
webpage19 <- read_html(url19)
description_data_html19 <- html_nodes(webpage19,'p')
description_data19 <- html_text(description_data_html19)
mydfm19 <- dfm(description_data19, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

#
url20 <- 'http://www.migrationpolicy.org/article/straddling-two-worlds-highly-skilled-migrants-senegambia-and-switzerland'
webpage20 <- read_html(url20)
description_data_html20 <- html_nodes(webpage20,'p')
description_data20 <- html_text(description_data_html20)
mydfm20 <- dfm(description_data20, remove = stopwords("english"), stem = TRUE, remove_punct = TRUE)

## Combining Dataset
MYDFM <- rbind(mydfm, mydfm1, mydfm2, mydfm3, mydfm4, mydfm5, mydfm6, mydfm7, mydfm8, mydfm9, mydfm10, mydfm11, mydfm12, mydfm13, mydfm14, mydfm15, mydfm16, mydfm17, mydfm18, mydfm19, mydfm20)

# Wordcloud
textplot_wordcloud(MYDFM, min.freq = 6,random.order = FALSE, rot.per = .25, colors = RColorBrewer::brewer.pal(8,"Dark2"))

## Sentiment Analysis
library(tidytext)
library(tidyverse)
library(stringr)
library(tidytext)

get_sentiments("afinn")
get_sentiments("bing")
get_sentiments("nrc")

titles <- c("description_data", "description_data1", "description_data2", "description_data3", "description_data4", "description_data5", "description_data6", "description_data7", "description_data8", "description_data9", "description_data10", "description_data11", "description_data12", "description_data13", "description_data14", "description_data15", "description_data16", "description_data17", "description_data18", "description_data19")

books <- list(description_data, description_data1, description_data2, description_data3, description_data4, description_data5,  description_data6, description_data7, description_data8, description_data9, description_data10, description_data11, description_data12, description_data13, description_data14, description_data15, description_data16, description_data17, description_data18,  description_data19)

series <- tibble()

for(i in seq_along(titles)) {

clean <- tibble(chapter = seq_along(books[[i]]), text = books[[i]]) %>%
unnest_tokens(word, text) %>%
mutate(book = titles[i]) %>%
select(book, everything())

series <- rbind(series, clean)
}

series$book <- factor(series$book, levels = rev(titles)) series %>%
right_join(get_sentiments("nrc")) %>%
filter(!is.na(sentiment)) %>%
count(sentiment, sort = TRUE)

library(ggplot2)

series %>%
group_by(book) %>%
mutate(word_count = 1:n(),
index = word_count %/% 100 + 1) %>%
inner_join(get_sentiments("bing")) %>%
count(book, index = index , sentiment) %>%
ungroup() %>%
spread(sentiment, n, fill = 0) %>%
mutate(sentiment = positive - negative,
book = factor(book, levels = titles)) %>%
ggplot(aes(index, sentiment, fill = book)) +
geom_bar(alpha = 0.5, stat = "identity", show.legend = FALSE) +
facet_wrap(~ book, ncol = 2, scales = "free_x")

# Comparing sentiment with different dictionaries (binds = 300 words)
afinn <- series %>%
group_by(book) %>%
mutate(word_count = 1:n(),
index = word_count %/% 300 + 1) %>%
inner_join(get_sentiments("afinn")) %>%
group_by(book, index) %>%
summarise(sentiment = sum(score)) %>%
mutate(method = "AFINN")

bing_and_nrc <- bind_rows(series %>%
group_by(book) %>%
mutate(word_count = 1:n(),
index = word_count %/% 300 + 1) %>%
inner_join(get_sentiments("bing")) %>%
mutate(method = "Bing"),
series %>%
group_by(book) %>%
mutate(word_count = 1:n(),
index = word_count %/% 300 + 1) %>%
inner_join(get_sentiments("nrc") %>%
filter(sentiment %in% c("positive", "negative"))) %>%
mutate(method = "NRC")) %>%
count(book, method, index = index , sentiment) %>%
ungroup() %>%
spread(sentiment, n, fill = 0) %>%
mutate(sentiment = positive - negative) %>%
select(book, index, method, sentiment)

bind_rows(afinn,
bing_and_nrc) %>%
ungroup() %>%
mutate(book = factor(book, levels = titles)) %>%
ggplot(aes(index, sentiment, fill = method)) +
geom_bar(alpha = 0.8, stat = "identity", show.legend = FALSE) +
facet_grid(book ~ method)

library(stm)
stm.content <- stm(MYDFM, K = 10)
labelTopics(stm.content)
topicCorr(stm.content)
plot(stm.content,type="summary")

library(plyr)
library(stringr)
library(e1071)

afinn_list <- read.delim(file='AFINN-111.txt', header=FALSE, stringsAsFactors=FALSE)
names(afinn_list) <- c('word', 'score')
afinn_list$word <- tolower(afinn_list$word)

posText <- read.delim(file='pos.txt', header=FALSE, stringsAsFactors=FALSE)
posText <- posText$V1
posText <- unlist(lapply(posText, function(x) { str_split(x, "\n") }))
negText <- read.delim(file='Neg.txt', header=FALSE, stringsAsFactors=FALSE)
negText <- negText$V1
negText <- unlist(lapply(negText, function(x) { str_split(x, "\n") }))

# Build tables of positive and negative sentences with scores
posResult <- as.data.frame(sentimentScore(posText, vNegTerms, negTerms, posTerms, vPosTerms))
negResult <- as.data.frame(sentimentScore(negText, vNegTerms, negTerms, posTerms, vPosTerms))
posResult <- cbind(posResult, 'positive')
colnames(posResult) <- c('sentence', 'vNeg', 'neg', 'pos', 'vPos', 'sentiment')
negResult <- cbind(negResult, 'negative')
colnames(negResult) <- c('sentence', 'vNeg', 'neg', 'pos', 'vPos', 'sentiment')

Labels:   train = [("restrict", "neg"), ("conflict", "neg"), ("ban", "neg"), ("deport", "neg"), ("barrier", "neg"),("agreement", "pos"),("assimilation", "pos"),("accept", "pos"),("connection", "pos")]

vNegTerms <- afinn_list$word[afinn_list$score==-5 | afinn_list$score==-4]
negTerms <- c(afinn_list$word[afinn_list$score==-3 | afinn_list$score==-2 | afinn_list$score==-1], "restrict", "conflict", "ban", "deport", "barrier")
posTerms <- c(afinn_list$word[afinn_list$score==3 | afinn_list$score==2 | afinn_list$score==1], "agreement", "assimilation", "accept", "connection")
vPosTerms <- c(afinn_list$word[afinn_list$score==5 | afinn_list$score==4])

# Functions to calculate numbers in each sentences.
SentimentScore <- function(sentences, vNegTerms, negTerms, posTerms, vPosTerms){
final_scores <- matrix('', 0, 5)
scores <- laply(sentences, function(sentence, vNegTerms, negTerms, posTerms, vPosTerms){
initial_sentence <- sentence
#Remove unnecessary characters and split up by word
sentence <- gsub('[[:punct:]]', '', sentence)
sentence <- gsub('[[:cntrl:]]', '', sentence)
sentence <- gsub('\\d+', '', sentence)
sentence <- tolower(sentence)
wordList <- str_split(sentence, '\\s+')
words <- unlist(wordList)

#Build vector with matches between sentence and each category
vPosMatches <- match(words, vPosTerms)
posMatches <- match(words, posTerms)
vNegMatches <- match(words, vNegTerms)
negMatches <- match(words, negTerms)

#Sum up number of words in each category
vPosMatches <- sum(!is.na(vPosMatches))
posMatches <- sum(!is.na(posMatches))
vNegMatches <- sum(!is.na(vNegMatches))
negMatches <- sum(!is.na(negMatches))
score <- c(vNegMatches, negMatches, posMatches, vPosMatches)

#Add row to scores table
newrow <- c(initial_sentence, score)
final_scores <- rbind(final_scores, newrow)
return(final_scores)
}, vNegTerms, negTerms, posTerms, vPosTerms)
return(scores)
}

#Build tables of positive and negative sentences with scores
posResult <- as.data.frame(sentimentScore(posText, vNegTerms, negTerms, posTerms, vPosTerms))
negResult <- as.data.frame(sentimentScore(negText, vNegTerms, negTerms, posTerms, vPosTerms))
posResult <- cbind(posResult, 'positive')
colnames(posResult) <- c('sentence', 'vNeg', 'neg', 'pos', 'vPos', 'sentiment')
negResult <- cbind(negResult, 'negative')
colnames(negResult) <- c('sentence', 'vNeg', 'neg', 'pos', 'vPos', 'sentiment')

#Combine the positive and negative tables
results <- rbind(posResult, negResult)

#Run the naive bayes algorithm using all four categories
classifier <- naiveBayes(results[,2:5], results[,6])

matrix <- table(predict(classifier, results), results[,6], dnn=list('predicted','actual'))
matrix

#Run a binomial test for confidence interval of results
binom.test(confTable[1,1] + confTable[2,2], nrow(results), p=0.5)

