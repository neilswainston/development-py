import feedparser as f
import markovify as m

_T = [[e.title, e.summary]
      for e in f.parse('feed://www.dailymail.co.uk/articles.rss').entries]

for _ in range(100):
    print(m.Text([j for sub in _T for j in sub]).make_short_sentence(140))
