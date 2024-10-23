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


<style>

    /* COLLAPSIBLE */
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

    /* Table */
    table, th, td {
    color: var(--global-text-color) !important; 

    table {
      margin: 0 auto; /* Center the table */
      width: auto; /* Let the table size according to content */
      text-align: left;
      border-spacing: 10px;
    }

    /* Adjust for smaller screens (two columns layout) */
    @media only screen and (max-width: 768px) {
      table {
        width: 100%; /* Ensure the table takes full width on mobile */
      }

      tr {
        display: flex;
        flex-wrap: wrap;
        justify-content: space-between;
      }

      td, th {
        display: inline-block;
        width: 45%; /* Two columns, 45% each */
        margin-bottom: 10px;
      }

      /* Give the last td more width */
      td:last-child {
        width: 100%; /* Allow the last column to take full width */
      }
    }
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

I am a PhD Candidate at Universitat Autònoma de Barcelona (UAB) and Barcelona School of Economics (BSE) in the [International Doctorate of Economic Analysis (IDEA) program](https://www.uabidea.eu). My supervisors are [Joan Llull](https://joanllull.github.io) and [Hanna Wang](https://sites.google.com/view/hannawang/).

**Research Interests**: Applied Microeconomics; Labor; Migration; Family Economics; Economic Geography

I'm on the Econ Job Market 2024. 

**Job Market Paper**:

[<u>"Parental Rural-Urban Migration and Child Education"</u>](https://conghanzheng.github.io/assets/pdf/Conghan_JMPaper2024.pdf)

<!-- 
Abstract [<button data-toggle="collapse" data-target="#aio" class="collapsible custom-collapsible"></button>]
<div id="aio" class="collapse" style="padding-left: 10px">
    Mobility frictions affect how much and where workers migrate and whether they bring children. This leads migrants to crowd into the most congested areas.
</div>
-->

**Job Market References**:

<table style="border-spacing:10px;">
  <tr>
    <td>
      <a href="https://joanllull.github.io">Joan Llull</a><br>
      <a href="mailto:joan.llull@iae.csic.es">joan.llull@iae.csic.es</a>
    </td>
    <td>
      <a href="https://sites.google.com/view/hannawang/">Hanna Wang</a><br>
      <a href="mailto:hanna.wang@uab.cat">hanna.wang@uab.cat</a>
    </td>
    <td>
      <a href="https://sites.google.com/view/adaferrer-i-carbonell">Ada Ferrer-i-Carbonell</a><br>
      <a href="mailto:ada.ferrer@iae.csic.es">ada.ferrer@iae.csic.es</a>
    </td>
  </tr>
</table>

<span style="display: block; margin-top: 50px;"></span>