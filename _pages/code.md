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

<span style="display: block; margin-top: 20px;"></span>

Below are some of the scripts I’ve developed for collecting online labor data and playing around with data. For course-related codes, please visit [my teaching page](https://conghanzheng.github.io/teaching/).

<span style="display: block; margin-top: 20px;"></span>

<ul>
    <li><b>Automated Data Collection</b> <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/5d2547cf7402bba5b8c6fd3c8a1c3e16d61a7318/assets/R/eurostatB.R">[R code]</a></li>
    <li><b>Applied Microeconometrics</b></li> 
        <ul>
            <li><b>Demand Estimation</b></li>
                <ul>
                    <li><b>Berry, Levinsohn & Pakes, 1995 (BLP); Nevo, 2000</b> <a href="https://conghanzheng.github.io/Applied_IO_TA/BLP.html">[Python Notebook]</a></li>
                    <li><b>Almost Ideal Demand System (AIDS)</b> <a href="https://conghanzheng.github.io/Applied_IO_TA/AIDS.html">[R Markdown]</a></li>
                </ul>
        </ul>
    <li><b>Other Econometric Topics</b></li> 
        <ul>
            <li><b>Kronecker Product Structure (KPS) Covariance</b> [Python Code<button data-toggle="collapse" data-target="#kps" class="collapsible custom-collapsible"></button>]
            <div id="kps" class="collapse">
            <span style="display: block; margin-top: 10px;"></span>
                <a href="https://github.com/conghanzheng/KPS">
                <img src="https://github-readme-stats.vercel.app/api/pin/?username=conghanzheng&repo=KPS&theme=transparent" alt="Readme Card">
                </a>
            </div></li>
        </ul>
</ul>