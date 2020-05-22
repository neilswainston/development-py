import feedparser as f
import markovify as m

_TEXT = [[e.title, e.summary]
         for e in f.parse('feed://www.dailymail.co.uk/articles.rss').entries]

for _ in range(100):
    print(m.Text([j for sub in _TEXT for j in sub]).make_short_sentence(140))
