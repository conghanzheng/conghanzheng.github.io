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

Below are some of the scripts I’ve developed for collecting online labor data and playing around with data. For other course-related codes, please visit [my teaching page](https://conghanzheng.github.io/teaching/).

<span style="display: block; margin-top: 20px;"></span>

<ul>
    <li><a>Automated Data Collection</a> <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/R/eurostatB.R">[R code]</a></li>
    <li><a>Applied Microeconometrics</a></li> 
        <ul>
            <li><a>Panel Data</a>: <a>Static and Dynamic Models</a> <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/R/PanelData.R">[R Code]</a> </li>
            <li><a>Discrete Choice</a>: <a>Binary and Multinomial Models</a> <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/R/DiscreteChoice.R">[R Code]</a> </li>
            <li><a>Censoring, Truncation, and Selection</a> <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/R/Selection.R">[R Code]</a> </li>
            <li><a>Demand Estimation</a></li>
                <ul>
                    <li><a>Berry, Levinsohn & Pakes, 1995 (BLP); Nevo, 2000</a> <a href="https://conghanzheng.github.io/Applied_IO_TA/BLP.html">[Python Notebook]</a></li>
                    <li><a>Almost Ideal Demand System (AIDS)</a> <a href="https://conghanzheng.github.io/Applied_IO_TA/AIDS.html">[R Markdown]</a></li>
                </ul>
        </ul>
    <li><a>Other Econometric Topics</a></li> 
        <ul>
            <li><a>Kronecker Product Structure (KPS) Covariance</a> <a href="https://github.com/conghanzheng/KPS">[Python Code]</a></li>
        </ul>
</ul>

<span style="display: block; margin-top: 50px;"></span>