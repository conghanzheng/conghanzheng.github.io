---
layout: page
permalink: /teaching/
title: Teaching
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

<span style="display: block; margin-top: 10px;"></span>

### Teaching Assistant
(Graduate Level Courses at IDEA Program, UAB)

<span style="display: block; margin-top: 10px;"></span>

**Microeconometrics, Part I** (Fall 2021, Fall 2022, Lecturer: Hanna Wang; Fall 2024, Lecturer: Joan Llull) [<button data-toggle="collapse" data-target="#mc1" class="collapsible custom-collapsible"></button>]
<div id="mc1" class="collapse">
   <ul>
    <li><a href="https://conghanzheng.github.io/assets/TA/MicroecomometricsI_2024/Syllabus.pdf">[Syllabus]</a> and <a href="https://joanllull.github.io/teaching.htm">Course Materials</a></li>
    <li><b>Panel Data</b></li>
    <ul>
        <li>Problem Set 1 (Due Sep 27) <a href="https://conghanzheng.github.io/assets/TA/MicroecomometricsI_2024/PS1.pdf">[PS1.pdf]</a> <a href="https://conghanzheng.github.io/assets/TA/MicroecomometricsI_2024/PS1_1.dta">[PS1_1.dta]</a> <a href="https://conghanzheng.github.io/assets/TA/MicroecomometricsI_2024/PS1_2.dta">[PS1_2.dta]</a> </li>
        <li>TA Session 1 <a href="https://conghanzheng.github.io/assets/TA/MicroecomometricsI_2024/TA1.pdf">[TA1.pdf]</a> <a href="https://conghanzheng.github.io/assets/TA/MicroecomometricsI_2024/TA1.do">[TA1.do]</a> <a href="https://conghanzheng.github.io/assets/TA/MicroecomometricsI_2024/TA1.dta">[TA1.dta]</a> </li>
    </ul>
    <li><b>Discrete Choice</b></li>
    <ul>
        <li>Problem Set 2 (Due Oct 04)</li>
    </ul>
    <li><b>Censoring, Truncation, and Selection</b></li>
    <ul>
        <li>Problem Set 3 (Due Oct 16)</li>
    </ul>
    <li><b>Duration Analysis</b></li>
    <ul>
        <li>Problem Set 4 (Due Nov 01)</li>
    </ul>
  </ul>
</div>

**Microeconometrics, Part II** (Fall 2024, Lecturer: Joan Llull) 

**Econometrics II** (Spring 2024, Lecturer: [Michael Creel](http://pareto.uab.es/mcreel/)) [<button data-toggle="collapse" data-target="#ec2" class="collapsible custom-collapsible"></button>]
<div id="ec2" class="collapse">
  <ul>
    <li><a href="https://conghanzheng.github.io/assets/TA/EcomometricsII_2024/Syllabus.pdf">[Syllabus]</a> and <a href="https://github.com/mcreel/Econometrics">Course Page</a></li>
    <li> Numerical Optimization <a href="https://conghanzheng.github.io/assets/TA/EcomometricsII_2024/TA1.m">[TA1.m]</a></li>
    <li> MLE <a href="https://conghanzheng.github.io/assets/TA/EcomometricsII_2024/TA2.m">[TA2.m]</a> <a href="https://conghanzheng.github.io/assets/TA/EcomometricsII_2024/TA2.ipynb">[TA2.ipynb]</a></li>
    <li> GMM <a href="https://conghanzheng.github.io/assets/TA/EcomometricsII_2024/TA3.m">[TA3.m]</a></li>
  </ul>
</div>

**Development Economics** (Spring 2024; Lecturer: [Laurence Go](https://www.laurencego.com)) [<button data-toggle="collapse" data-target="#dev" class="collapsible custom-collapsible"></button>]
<div id="dev" class="collapse">
  <span>&nbsp;&nbsp;&nbsp;&nbsp;</span> <a href="https://conghanzheng.github.io/assets/TA/Dev_Econ_2024/Syllabus.pdf"><b> [Syllabus] </b></a> <br>
  <span>&nbsp;&nbsp;&nbsp;&nbsp;</span>
</div>

**Applied Industrial Organization** (Spring 2023, Lecturer: [Susanna Esteban](https://www.cemfi.es/people/faculty/profile.asp?u=esteban)) [<button data-toggle="collapse" data-target="#aio" class="collapsible custom-collapsible"></button>]
<div id="aio" class="collapse">
  <span>&nbsp;&nbsp;&nbsp;&nbsp;</span> <a href="https://conghanzheng.github.io/assets/TA/Applied_IO_2023/Syllabus.pdf"><b> [Syllabus] </b></a> <br>
  <span style="display: block; margin-top: 10px;"></span>
  <span>&nbsp;&nbsp;&nbsp;&nbsp;</span> <a href="https://github.com/conghanzheng/Applied_IO_TA">
        <img src="https://github-readme-stats.vercel.app/api/pin/?username=conghanzheng&repo=Applied_IO_TA&theme=transparent" alt="Readme Card">
  </a>
</div>