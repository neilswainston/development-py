import feedparser
import markovify

d = feedparser.parse('feed://www.dailymail.co.uk/articles.rss')
text = []

for post in d.entries:
    text.append(post.title)
    text.append(post.summary)

text_model = markovify.Text(text)

for i in range(100):
    print(text_model.make_short_sentence(140))
