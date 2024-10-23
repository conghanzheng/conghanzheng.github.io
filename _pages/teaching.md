---
layout: page
permalink: /teaching/
title: Teaching
description: <a>*Teaching evaluations available upon request.</a>
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

<br>

#### Teaching Assistant
(Graduate Level Courses at IDEA Program, UAB)

<span style="display: block; margin-top: 20px;"></span>

Microeconometrics, Part I (Fall 2021, Fall 2022, Lecturer: Hanna Wang; Fall 2024, Lecturer: Joan Llull) [<button data-toggle="collapse" data-target="#mc1" class="collapsible custom-collapsible"></button>]
<div id="mc1" class="collapse">
    <a>*Please download the problem set pdf, solution pdf and codes before <u>December 13, 2024</u>. After that, all problem sets files will be deleted from here.</a>
    <span style="display: block; margin-top: 10px;"></span>
    <a>*For the Stata codes, I adapt the scripts contributed by previous TAs (Manuel V. Montesinos, Chanon Techasunthornwat, and Bartosz Korchowiec). I also refer to the book Cameron and Trivedi (2022): "Microeconometrics using stata". All errors are mine.</a>
    <span style="display: block; margin-top: 10px;"></span>
   <ul>
    <li><a href="https://conghanzheng.github.io/assets/TA/MicroeconometricsI_2024/Syllabus.pdf">[Syllabus]</a> and <a href="https://joanllull.github.io/teaching.htm">Course Materials</a></li>
    <li><a>Panel Data</a></li>
    <ul>
        <li>Problem Set 1 (Due Sep 27) </li>
            <ul>
                <li><a href="https://conghanzheng.github.io/assets/TA/MicroeconometricsI_2024/PS1.pdf">[PS1.pdf]</a>, Data: <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/PS1_1.dta">[PS1_1.dta]</a> <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/PS1_1.dta">[PS1_2.dta]</a> </li>
                <li>Solution: <a href="https://conghanzheng.github.io/assets/TA/MicroeconometricsI_2024/PS1_solution.pdf">[PS1_solution.pdf]</a>, Code: <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/R/PanelData.R">[PS1.R]</a>, <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/PS1.do">[PS1.do]</a>, <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/PS1_ABover.m">[ABover.m]</a> </li>
            </ul>
        <li>TA Session (Sep 20)</li>
            <ul>
                <li>Slides: <a href="https://conghanzheng.github.io/assets/TA/MicroeconometricsI_2024/TA1.pdf">[TA1.pdf]</a></li>
                <li>Code: <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/TA1.do">[TA1.do]</a>, Data: <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/TA1.dta">[TA1.dta]</a></li>
            </ul>
    </ul>
    <li><a>Discrete Choice</a></li>
    <ul>
        <li>Problem Set 2 (Due Oct 04)</li>
            <ul>
                <li><a href="https://conghanzheng.github.io/assets/TA/MicroeconometricsI_2024/PS2.pdf">[PS2.pdf]</a>, Data: <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/PS2.dta">[PS2.dta]</a></li>
                <li>Solution: <a href="https://conghanzheng.github.io/assets/TA/MicroeconometricsI_2024/PS2_solution.pdf">[PS2_solution.pdf]</a>, Code: <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/R/DiscreteChoice.R">[PS2.R]</a>, <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/PS2.do">[PS2.do]</a></li>
            </ul>
        <li>TA Session (Sep 27)</li>
            <ul>
                <li>Slides:  <a href="https://conghanzheng.github.io/assets/TA/MicroeconometricsI_2024/TA2.pdf">[TA2.pdf]</a></li>
                <li>Code: <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/TA2.do">[TA2.do]</a>, Data: <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/TA2_1.dta">[TA2_1.dta]</a> <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/TA2_2.dta">[TA2_2.dta]</a> <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/TA2_3.dta">[TA2_3.dta]</a></li>
            </ul>
    </ul>
    <li><a>Censoring, Truncation, and Selection</a></li>
    <ul>
        <li>Problem Set 3 (Due Oct 16)</li>
            <ul>
                <li><a href="https://conghanzheng.github.io/assets/TA/MicroeconometricsI_2024/PS3.pdf">[PS3.pdf]</a>, Data: <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/PS3.dta">[PS3.dta]</a>, <i>Feedback on your solutions will be sent out during the week of Oct 21.</i></li>
                <li>Solution: <a href="https://conghanzheng.github.io/assets/TA/MicroeconometricsI_2024/PS3_solution.pdf">[PS3_solution.pdf]</a>, Code: <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/R/Selection.R">[PS3.R]</a>, <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/PS3.do">[PS3.do]</a></li>
            </ul>
        <li>TA Session (Oct 11)</li>
        <ul>
            <li>Slides:  <a href="https://conghanzheng.github.io/assets/TA/MicroeconometricsI_2024/TA3.pdf">[TA3.pdf]</a></li>
            <li>Code: <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/TA3.do">[TA3.do]</a>, Data: <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/TA3.dta">[TA3.dta]</a></li>
        </ul>
    </ul>
    <li><a>Duration Analysis</a></li>
    <ul>
        <li>Problem Set 4 (Due Nov 01)</li>
            <ul>
                <li><a href="https://conghanzheng.github.io/assets/TA/MicroeconometricsI_2024/PS4.pdf">[PS4.pdf]</a>, Data: <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/PS4_1.dta">[PS4_1.dta]</a> <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/PS4_2.dta">[PS4_2.dta]</a>, <i>Feedback on your solutions will be sent out around Nov 05</i></li>
                <li>Solution: <i>will be posted on Nov 02</i></li>
            </ul>
        <li>TA Session (Oct 16):</li>
            <ul>
                <li>Slides:  <a href="https://conghanzheng.github.io/assets/TA/MicroeconometricsI_2024/TA4.pdf">[TA4.pdf]</a></li>
                <li>Code: <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/TA4.do">[TA4.do]</a>, Data: <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/MicroeconometricsI_2024/TA4.dta">[TA4.dta]</a></li>
            </ul>
    </ul>
  </ul>
</div>

Microeconometrics, Part II (Fall 2024, Lecturer: Joan Llull) [<button data-toggle="collapse" data-target="#mc2" class="collapsible custom-collapsible"></button>]
<div id="mc2" class="collapse">
    <i>TA Sessions</i> : <a>We will have a total of four sessions, each session will cover one of the following four topics. For each topic, I will post tutorial TA codes and slides. </a>
    <span style="display: block; margin-top: 10px;"></span>
    <i>Problem Sets</i> : <a>The solutions to the first two problem sets will be discussed in class. For all problem sets, I will post solution codes. Feedback on your solution will be provided by marking and commenting on your solution PDF, which will be sent to you within one week after it is due.</a>
     <span style="display: block; margin-top: 10px;"></span>
    <a>*Please download the problem set pdf and solution code scripts before <u>December 13, 2024</u>. After that, all problem sets files will be deleted from here.</a>
    <span style="display: block; margin-top: 10px;"></span>
    <ul>
        <li><a href="https://conghanzheng.github.io/assets/TA/MicroeconometricsII_2024/Syllabus.pdf">[Syllabus]</a> and <a href="https://joanllull.github.io/teaching.htm">Course Materials</a></li>
        <li><a>Treatment Effects: RCT, Matching, and IV</a></li>
            <ul>
                <li>Problem Set 1 (Due Nov 15) </li>
                <li>TA Session (Nov 08)</li>
            </ul>
        <li><a>Treatment Effects: RD and DiD</a></li>
            <ul>
                <li>Problem Set 2 (Due Nov 22) </li>
                <li>TA Session (Nov 15)</li>
            </ul>
        <li><a>Dynamic Discrete Choice: Full Solution</a></li>
            <ul>
                <li>Problem Set 3 </li>
                <li>TA Session (Nov 22)</li>
            </ul>
        <li><a>Dynamic Discrete Choice: CCP</a></li>
            <ul>
                <li>Problem Set 4 </li>
                <li>TA Session (Nov 29)</li>
            </ul>
    </ul>
</div>

Econometrics II (Spring 2024, Lecturer: [Michael Creel](http://pareto.uab.es/mcreel/)) [<button data-toggle="collapse" data-target="#ec2" class="collapsible custom-collapsible"></button>]
<div id="ec2" class="collapse">
  <ul>
    <li><a href="https://conghanzheng.github.io/assets/TA/EconometricsII_2024/Syllabus.pdf">[Syllabus]</a> and <a href="https://github.com/mcreel/Econometrics">Course Page</a></li>
    <li> Numerical Optimization <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/EconometricsII_2024/TA1.m">[TA1.m]</a></li>
    <li> MLE <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/EconometricsII_2024/TA2.m">[TA2.m]</a> <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/EconometricsII_2024/TA2.ipynb">[TA2.ipynb]</a></li>
    <li> GMM <a href="https://github.com/conghanzheng/conghanzheng.github.io/blob/master/assets/TA/EconometricsII_2024/TA3.m">[TA3.m]</a></li>
  </ul>
</div>

Development Economics (Spring 2024; Lecturer: [Laurence Go](https://www.laurencego.com)) [<button data-toggle="collapse" data-target="#dev" class="collapsible custom-collapsible"></button>]
<div id="dev" class="collapse">
  <ul>
    <li><a href="https://conghanzheng.github.io/assets/TA/Dev_Econ_2024/Syllabus.pdf">[Syllabus]</a></li>
  </ul>
</div>

Applied Industrial Organization (Spring 2023, Lecturer: [Susanna Esteban](https://www.cemfi.es/people/faculty/profile.asp?u=esteban)) [<button data-toggle="collapse" data-target="#aio" class="collapsible custom-collapsible"></button>]
<div id="aio" class="collapse">
    <ul>
        <li><a href="https://conghanzheng.github.io/assets/TA/Applied_IO_2023/Syllabus.pdf">[Syllabus]</a></li>
    </ul>
  <span>&nbsp;&nbsp;&nbsp;&nbsp;</span> <a href="https://github.com/conghanzheng/Applied_IO_TA">
        <img src="https://github-readme-stats.vercel.app/api/pin/?username=conghanzheng&repo=Applied_IO_TA&theme=transparent" alt="Readme Card">
  </a>
</div>

<span style="display: block; margin-top: 50px;"></span>