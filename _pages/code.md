---
layout: page
permalink: /code/
title: Code
description: 
nav: true
nav_order: 3
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

<br>

Below are some of the scripts I’ve developed for collecting online labor data and playing around with data. For course-related codes, please visit [my teaching page](https://conghanzheng.github.io/teaching/).

<span style="display: block; margin-top: 20px;"></span>

<ul>
    <li><a>Automated Data Collection</a> <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/5d2547cf7402bba5b8c6fd3c8a1c3e16d61a7318/assets/R/eurostatB.R">[R code]</a></li>
    <li><a>Applied Microeconometrics</a></li> 
        <ul>
            <li><a>Panel Data</a>: <a>Static and Dynamic</a> <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/fe155a9a3e160738cd6030dcb64d419c27c463a3/assets/TA/MicroeconometricsI_2024/PS1.R">[R Code]</a> </li>
            <li><a>Demand Estimation</a></li>
                <ul>
                    <li><a>Berry, Levinsohn & Pakes, 1995 (BLP); Nevo, 2000</a> <a href="https://conghanzheng.github.io/Applied_IO_TA/BLP.html">[Python Notebook]</a></li>
                    <li><a>Almost Ideal Demand System (AIDS)</a> <a href="https://conghanzheng.github.io/Applied_IO_TA/AIDS.html">[R Markdown]</a></li>
                </ul>
        </ul>
    <li><a>Other Econometric Topics</a></li> 
        <ul>
            <li><a>Kronecker Product Structure (KPS) Covariance</a> [Python Code<button data-toggle="collapse" data-target="#kps" class="collapsible custom-collapsible"></button>]
            <div id="kps" class="collapse">
            <span style="display: block; margin-top: 10px;"></span>
                <a href="https://github.com/conghanzheng/KPS">
                <img src="https://github-readme-stats.vercel.app/api/pin/?username=conghanzheng&repo=KPS&theme=transparent" alt="Readme Card">
                </a>
            </div></li>
        </ul>
</ul>