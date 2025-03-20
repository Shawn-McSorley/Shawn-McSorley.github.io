---
layout: default
title: Home
---
# My notes on software defined radio

{% if site.tags.SDR %}
  <ul>
    {% for post in site.tags.SDR %}
      <li><a href="{{ post.url }}">{{ post.title }}</a></li>
    {% endfor %}
  </ul>
{% else %}
  <p>No posts found under 'SDR'.</p>
{% endif %}
# Latest Posts
<ul>
  {% for post in site.posts limit: 3 %}
    <li>
      <a href="{{ post.url }}">{{ post.title }}</a> - {{ post.date | date: "%B %d, %Y" }}
    </li>
  {% endfor %}
</ul>


