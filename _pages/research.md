---
layout: page
permalink: /research/
title: Research
description: 
nav: true
nav_order: 1
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

<span style="display: block; margin-top: 30px;"></span>

#### Job Market Paper

<span style="display: block; margin-top: 20px;"></span>

"Parental Rural-Urban Migration and Child Education", 2024

<span style="display: block; margin-top: 30px;"></span>

#### Working Paper

<span style="display: block; margin-top: 20px;"></span>

[<u>"Migration Outcome Gap: The Cost of Leaving Children Behind"</u>](https://conghanzheng.github.io/assets/pdf/ENTER_JBR_paper.pdf), 2023 <br>
<span>&nbsp;</span> Abstract [<button data-toggle="collapse" data-target="#enterjbr2023" class="collapsible custom-collapsible"></button>]
<div id="enterjbr2023" class="collapse">
    <p>
        Parents invest money and parental care in their children. Given the constraints on migration, migrant parents may trade one for the other. In China, in addition to financial constraints, there are also institutional mobility restrictions on rural-urban migration, which limit the ability of migrants to claim urban citizenship and thus pose the problem of multiple selectivity. Using a simultaneous equations model to deal with endogeneity, the results highlight the motive of households to migrate for better educational opportunities for children and predict a gain from rural-urban migration in children's educational outcomes. And a significant part of the income gains from migration is invested in the child's human capital accumulation. The results also suggest that <i>Hukou</i> restrictions in China don't prevent migrant families from going to urban areas, but only limit their access to local social services and lead to a loss of welfare.
    </p>
</div>

<span style="display: block; margin-top: 30px;"></span>

#### Work-in-Progress

<span style="display: block; margin-top: 20px;"></span>

"Family Migration and the Intergenerational Contract"

<span style="display: block; margin-top: 30px;"></span>

#### Publication

<span style="display: block; margin-top: 20px;"></span>

[<u>"Forecasting Bilateral Refugee Flows with High-dimensional Data and Machine Learning Techniques"</u>](https://doi.org/10.1093/jeg/lbae023) (with [Konstantin Boss](https://sites.google.com/view/konstantinboss/home) and [Tobias Heidland](https://www.tobias-heidland.de) and [Andre Groeger](https://sites.google.com/site/andregroeger/) and Finja Krüger) <br>
*Journal of Economic Geography*, August 2024: lbae023.<br>
<span>&nbsp;</span> Abstract [<button data-toggle="collapse" data-target="#jeg2024" class="collapsible custom-collapsible"></button>]
<div id="jeg2024" class="collapse">
    <p>
        We develop monthly asylum seeker flow forecasting models for 157 origin countries to the EU27, using machine learning and high-dimensional data, including digital trace data from Google Trends. Comparing different models and forecasting horizons and validating out-of-sample, we find that an ensemble forecast combining Random Forest and Extreme Gradient Boosting algorithms outperforms the random walk over horizons between 3 and 12 months. For large corridors, this holds in a parsimonious model exclusively based on Google Trends variables, which has the advantage of near real-time availability. We provide practical recommendations how our approach can enable ahead-of-period asylum seeker flow forecasting applications.
    </p>
</div>

<br>