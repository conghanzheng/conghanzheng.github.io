---
layout: page
permalink: /teaching/
title: Teaching
description: 
nav: true
nav_order: 2
---

---
<br>

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
        content: "\2212";
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

#### Teaching Assistant, PhD Level Courses at IDEA Program, UAB

- Econometrics II (Spring 2024; Lecturer: [Michael Creel](http://pareto.uab.es/mcreel/)) [<button data-toggle="collapse" data-target="#ec2" class="collapsible"></button>]
<div id="ec2" class="collapse">
  <ul>
    <li><a href="https://drive.google.com/file/d/1GSNJRoYPvwdxkAYcsPvlopU8SRnx1VmM/view?usp=share_link">Syllabus</a> and <a href="https://github.com/mcreel/Econometrics">Course Page</a></li>
    <li>TA Session: Numerical Optimization <a href="https://drive.google.com/file/d/1QFRB7ZZL_lyMpmfvHX9CT163sfn2Dg5x/view?usp=share_link">[TA1.m]</a></li>
    <li>TA Session: MLE <a href="https://drive.google.com/file/d/1t7izaSqsrUKdleLYzHuhQh7LAU9qx0kF/view?usp=share_link">[TA2.m]</a> <a href="https://drive.google.com/file/d/1jMbmNokGRmbOs-8ZcRjyQRDYkZeUPD2y/view?usp=share_link">[TA2.ipynb]</a></li>
    <li>TA Session: GMM <a href="https://drive.google.com/file/d/1Bjhtpc_YfR9yMyIQVeM6VMPPiXTOhDn_/view?usp=share_link">[TA3.m]</a></li>
  </ul>
</div>

- [Microeconometrics I](https://drive.google.com/file/d/121f153kIO8coppfRyuKmPNtTCI1Kf7sP/view?usp=share_link) (Fall 2021, Fall 2022, Lecturer: Hanna Wang; Fall 2024, Lecturer: Joan Llull) [<button data-toggle="collapse" data-target="#mc1" class="collapsible"></button>]
<div id="mc1" class="collapse">
  <a href="https://github.com/conghanzheng/MicroeconometricsI_TA">
        <img src="https://github-readme-stats.vercel.app/api/pin/?username=conghanzheng&repo=MicroeconometricsI_TA&theme=transparent" alt="Readme Card">
  </a>
</div>

- [Microeconometrics II](https://drive.google.com/file/d/1D-mKQAbUQeMbdeWGhU87N7VMlOTwSIl1/view?usp=share_link) (Fall 2024; Lecturer: Joan Llull) 

- Development Economics (Spring 2024; Lecturer: [Laurence Go](https://www.laurencego.com))

- [Applied Industrial Organization](https://drive.google.com/file/d/1_Mo3X_meH9c37PkeOMv6Y6OaQCaE9kO7/view?usp=share_link) (Spring 2023; Lecturer: [Susanna Esteban](https://www.cemfi.es/people/faculty/profile.asp?u=esteban)) [<button data-toggle="collapse" data-target="#aio" class="collapsible"></button>]
<div id="aio" class="collapse">
  <a href="https://github.com/conghanzheng/Applied_IO_TA">
        <img src="https://github-readme-stats.vercel.app/api/pin/?username=conghanzheng&repo=Applied_IO_TA&theme=transparent" alt="Readme Card">
  </a>
</div>