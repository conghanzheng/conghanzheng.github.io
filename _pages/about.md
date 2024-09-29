---
layout: about
title: Home
permalink: /
subtitle: > # <p> (pronounced: <em>tsong-hahn</em>) </p>
  <a>PhD Candidate, UAB and BSE</a> <br>
  <span>&#9993;</span> <a href="mailto:zhengconghan@outlook.com">zhengconghan@outlook.com</a>

profile:
  align: right
  image: conghan2022.jpg
  image_circular: false # crops the image to make it circular
  more_info:

news: false # includes a list of news items
selected_papers: false # includes a list of papers marked as "selected={true}"
social: false # includes social icons at the bottom of the page
---


<!--
<i class="ai ai-orcid" style="color:#A6CE39;"></i> [ORCiD](https://orcid.org/0000-0003-0158-5111)
-->

<!-- COLLAPSIBLE-->
<style>
    .collapsible {
        background-color: transparent;
        cursor: pointer;
        border: none;
    }

    .collapsible:after {
        color: #6c0a31;
        content: '\002B'; /* Default plus sign */
    }

    .active:after {
        color: #6c0a31;
        content: ""; /* Minus sign for active state */
    }

    .custom-collapsible:after {
        content: '\002B'; /* Default plus sign for custom collapsible */
    }

    .custom-collapsible.active:after {
        content: "\2212"; /* Minus sign for active custom collapsible */
    }

    .content {
        max-height: 0;
        overflow: hidden;
        transition: max-height 0.2s ease-out;
    }
</style>

<script>
    document.addEventListener("DOMContentLoaded", function() {
        var coll = document.getElementsByClassName("collapsible");
        for (var i = 0; i < coll.length; i++) {
            coll[i].addEventListener("click", function() {
                this.classList.toggle("active");
                var content = this.nextElementSibling;
                if (content.style.maxHeight) {
                    content.style.maxHeight = null;
                } else {
                    content.style.maxHeight = content.scrollHeight + "px";
                }
            });
        }
    });
</script>

<span style="display: block; margin-top: 20px;"></span>

I am a PhD Candidate at Universitat Autònoma de Barcelona (UAB) and Barcelona School of Economics (BSE) in the [International Doctorate of Economic Analysis (IDEA) program](https://www.uabidea.eu). My supervisors are [Joan Llull](https://joanllull.github.io) and [Hanna Wang](https://sites.google.com/view/hannawang/).

**Research Interests**: Applied Microeconomics; Labor; Migration; Family Economics; Economic Geography

I'm on the Econ Job Market 2024. My placement director is [Inés Macho-Stadler](https://www.inesmachostadler.com).

<span style="display: block; margin-top: 20px;"></span>