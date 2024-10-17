---
layout: about
title: Home
permalink: /
subtitle: > # <p> (pronounced: <em>tsong-hahn</em>) </p>
  <a>PhD Candidate, UAB and BSE</a> <br>
  <span>&#9993;</span> <a href="mailto:zhengconghan@outlook.com">zhengconghan@outlook.com</a>

profile:
  align: right
  image: conghan2024.JPG
  image_circular: false # crops the image to make it circular
  more_info:

news: false # includes a list of news items
selected_papers: false # includes a list of papers marked as "selected={true}"
social: false # includes social icons at the bottom of the page
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

<span style="display: block; margin-top: 30px;"></span>

<span style="margin-top: 30px;"></span>

I am a PhD Candidate at Universitat Autònoma de Barcelona (UAB) and Barcelona School of Economics (BSE) in the [International Doctorate of Economic Analysis (IDEA) program](https://www.uabidea.eu). My supervisors are [Joan Llull](https://joanllull.github.io) and [Hanna Wang](https://sites.google.com/view/hannawang/).

**Research Interests**: Applied Microeconomics; Labor; Migration; Family Economics; Economic Geography

I'm on the Econ Job Market 2024. 

<!-- 

In my job market paper I study how rural-urban migrant workers decisions affect their children's educational outcomes under mobility constraints, and also I found that the migrants' location choice, their children's school location choice .

rural-urban migrants' decision on work location and their children's school location. 

[<u>"Parental Rural-Urban Migration and Child Education"</u>](https://conghanzheng.github.io/assets/pdf/Conghan_JMPaper2024.pdf), Job Market Paper, 2024 <br>
<i style="padding-left: 10px"> Abstract: </i> [<button data-toggle="collapse" data-target="#jmp2024" class="collapsible custom-collapsible"></button>]
<p style="padding-left: 10px" id="jmp2024" class="collapse"> 
    [Abstract Text Here]
</p>
-->

<span style="display: block; margin-top: 50px;"></span>

<!-- Add logos below 
<p align="center">
  <a href="https://www.uab.cat">
    <img src="/assets/img/UAB_Logo.png" width="30%" style="margin-right: 30px;" alt="UAB Logo">
  </a>
  <a href="https://www.uabidea.eu/home">
    <img src="/assets/img/IDEA_Logo.png" width="30%" style="margin-right: 30px;" alt="IDEA Logo">
  </a>
  <a href="https://bse.eu">
    <img src="/assets/img/BSE_Logo.jpg" width="30%" alt="BSE Logo">
  </a>
</p>
-->