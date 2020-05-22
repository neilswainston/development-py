import feedparser as f
import markovify as m

text = []

for ent in f.parse('feed://www.dailymail.co.uk/articles.rss').entries:
    text.extend([ent.title, ent.summary])

text_model = m.Text(text)

for _ in range(100):
    print(text_model.make_short_sentence(140))
