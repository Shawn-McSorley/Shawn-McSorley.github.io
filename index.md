---
layout: home
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