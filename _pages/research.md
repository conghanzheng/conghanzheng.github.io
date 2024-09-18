---
layout: page
permalink: /research/
title: Research
description:
nav: true
nav_order: 2
---

<!-- COLLAPSIBLE-->
<style>
    .collapsible {
        background-color: transparent;
        cursor: pointer;
        border: none;
    }

    .collapsible:after {
        color: #6c0a31;
        content: '\002B';
    }

    .active:after {
        color: #6c0a31;
        content: "";
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

---
<br>

### Job Market Paper


**Parental Rural-Urban Migration and Child Education**

<br>


### Publication


[**Forecasting Bilateral Refugee Flows with High-dimensional Data and Machine Learning Techniques**](https://doi.org/10.1093/jeg/lbae023) (with [Konstantin Boss](https://sites.google.com/view/konstantinboss/home) and [Tobias Heidland](https://www.tobias-heidland.de) and [Andre Groeger](https://sites.google.com/site/andregroeger/) and Finja Krüger), **Journal of Economic Geography**, August 2024

<span>&nbsp;</span> Abstract [<button data-toggle="collapse" data-target="#jeg2024" class="collapsible"></button>]
<div id="jeg2024" class="collapse">
    <p>
        We develop monthly asylum seeker flow forecasting models for 157 origin countries to the EU27, using machine learning and high-dimensional data, including digital trace data from Google Trends. Comparing different models and forecasting horizons and validating out-of-sample, we find that an ensemble forecast combining Random Forest and Extreme Gradient Boosting algorithms outperforms the random walk over horizons between 3 and 12 months. For large corridors, this holds in a parsimonious model exclusively based on Google Trends variables, which has the advantage of near real-time availability. We provide practical recommendations how our approach can enable ahead-of-period asylum seeker flow forecasting applications.
    </p>
</div>