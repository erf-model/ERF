.. _ERFKnowledgeAssistant:

ERF Knowledge Assistant
=======================

.. raw:: html

   <style>
   .erf-assistant-hero {
     margin: 1.25rem 0 1.5rem;
     padding: 1.25rem 1.4rem;
     border: 1px solid #d9e2ec;
     border-radius: 14px;
     background: linear-gradient(135deg, #f6fbff 0%, #eef6ff 48%, #f9fbff 100%);
     box-shadow: 0 10px 30px rgba(15, 23, 42, 0.08);
   }
   .erf-assistant-hero h2 {
     margin: 0 0 0.4rem;
     font-size: 1.45rem;
   }
   .erf-assistant-hero-context7 {
     background: linear-gradient(135deg, #f7fcf4 0%, #edf8ef 48%, #f8fcf7 100%);
     border-color: #d6ead8;
   }
   .erf-assistant-hero p {
     margin: 0.35rem 0 0;
     line-height: 1.5;
   }
   .erf-assistant-badge {
     display: inline-block;
     margin-bottom: 0.7rem;
     padding: 0.25rem 0.6rem;
     border-radius: 999px;
     background: #dbeafe;
     color: #1d4ed8;
     font-size: 0.82rem;
     font-weight: 700;
     letter-spacing: 0.02em;
     text-transform: uppercase;
   }
   .erf-assistant-actions {
     display: flex;
     flex-wrap: wrap;
     gap: 0.6rem;
     margin-top: 1rem;
   }
   .erf-assistant-button {
     display: inline-flex;
     align-items: center;
     justify-content: center;
     padding: 0.6rem 0.9rem;
     border: 1px solid #cbd5e1;
     border-radius: 999px;
     background: rgba(255, 255, 255, 0.82);
     color: #0f172a !important;
     font-size: 0.95rem;
     font-weight: 600;
     line-height: 1.2;
     text-decoration: none !important;
     transition: background-color 0.15s ease, border-color 0.15s ease, box-shadow 0.15s ease;
   }
   .erf-assistant-button:hover {
     background: #ffffff;
     border-color: #94a3b8;
     box-shadow: 0 6px 16px rgba(15, 23, 42, 0.08);
   }
   </style>

.. raw:: html

   <div class="erf-assistant-hero">
     <div class="erf-assistant-badge">Featured Assistant</div>
     <h2>ERF Specialized Assistant</h2>
     <p>
       Use this assistant to answer questions about the ERF code base, interpret
       documentation, and keep responses grounded in ERF-specific sources.
     </p>
     <p>
       If your browser does not support embeds for external services, use the
       button below.
     </p>
     <div class="erf-assistant-actions">
       <a class="erf-assistant-button" href="https://gemini.google.com/gem/1Kf8jWsF2tKLOyOzNDrjFkHsQNBD69mId?usp=sharing" target="_blank" rel="noopener noreferrer">
         Open the ERF Specialized Assistant
       </a>
     </div>
   </div>

.. raw:: html

   <div class="erf-assistant-hero erf-assistant-hero-context7">
     <div class="erf-assistant-badge">Context7 Chat</div>
     <h2>Context7 ERF Chatbot</h2>
     <p>
       Use the Context7 ERF chatbot for an MCP-backed ERF chat surface, or
       open the ERF and AMReX chatbots directly in Context7.
     </p>
     <div class="erf-assistant-actions">
       <a class="erf-assistant-button" href="https://context7.com/erf-model/erf?tab=chat" target="_blank" rel="noopener noreferrer">
         Open ERF Chatbot
       </a>
       <a class="erf-assistant-button" href="https://context7.com/amrex-codes/amrex?tab=chat" target="_blank" rel="noopener noreferrer">
         Open AMReX Chatbot
       </a>
     </div>
     <p>
       Context7 developer references:
       <a href="https://context7.com/docs/resources/developer" target="_blank" rel="noopener noreferrer">Developer docs</a>
       and
       <a href="https://context7.com/docs/llms.txt" target="_blank" rel="noopener noreferrer">documentation index</a>.
     </p>
     <script src="https://context7.com/widget.js" data-library="/erf-model/erf"></script>
   </div>

.. list-table::
   :widths: 30 18 18 18
   :header-rows: 1

   * - Assistant
     - Documented Here
     - MCP
     - ERF-specific
   * - `ERF Specialized Assistant <https://gemini.google.com/gem/1Kf8jWsF2tKLOyOzNDrjFkHsQNBD69mId?usp=sharing>`_
     - Yes
     - No
     - Yes
   * - | `Sphinx SingleHTML <https://gemini.google.com/gem/1wnuJ84MnXQrLnoZp3hvTCVsy_7tBmVck?usp=sharing>`_
       | Deep research (AMReX/ERF)
     - No
     - No
     - Yes
   * - | `Sphinx SingleHTML <https://gemini.google.com/gem/1-nyIxXiueWECypPG1QJrkTAuc2GTPUAb?usp=sharing>`_
       | Source/Exec
     - No
     - No
     - Yes
   * - `Context7 ERF chatbot <https://context7.com/erf-model/erf?tab=chat>`_
     - Yes
     - Yes
     - Yes
   * - `Context7 AMReX chatbot <https://context7.com/amrex-codes/amrex?tab=chat>`_
     - Yes
     - Yes
     - No

.. Previous assistant link retained for reference:
.. https://gemini.google.com/gem/1CYi43osCZtA-pqmuBOyw6AZQJpkQBHox?usp=sharing

ERF users can use a specialized assistant to answer questions about the code base,
help interpret documentation, and keep responses grounded in ERF-specific sources.
Use the button above to open the ready-to-use assistant.

This page collects the small, reusable assets that matter for ERF assistant setup.
It is intentionally separate from :doc:`AgenticWorkFlow`, which covers running
`AMReX-Agent <https://github.com/AMReX-Codes/amrex-agent>`_ with ERF.

How the assistant adapts
------------------------

The ERF assistant can adapt its level of detail to the user's background. At the
start of a session, it can ask whether you are new to HPC / AMReX / ERF,
familiar with HPC but newer to ERF, or an experienced ERF / AMReX user. That
answer shapes how much context, terminology, and step-by-step detail it provides.

It can also help author or modify ERF inputs files. In that mode, it should ask a
few clarifying questions first, then produce a commented template using only
documented parameter keys. Any unconfirmed key should be marked explicitly and
left for verification.

No advanced tools are enabled by default. If available in your assistant
environment, *Guided Learning* can convert the exchange into a structured,
step-by-step instructional workflow, while *Deep Research* can be used to
assemble a more evidence-heavy answer by gathering and synthesizing larger
amounts of documentation.

User expertise and assistant mode
---------------------------------

The assistant distinguishes between two separate forms of adaptation.

1. **User expertise leveling** controls the style of explanation. A novice user
   receives more definitions, more context, and more explicit step-by-step
   guidance. An experienced user receives a terser response with less framing.
   This concerns presentation, not the underlying evidence standard.
2. **Tool selection** controls how the assistant gathers or structures material.
   *Guided Learning* is suited to a structured instructional sequence, while
   *Deep Research* is suited to a broader evidence-gathering pass before the
   response is written. These tools are optional and not enabled by default.

The goal in both cases is to keep the assistant evidence-first: it should avoid
guessing ERF parameters, inventing undocumented keys, or claiming certainty
without supporting inputs.

What to use
-----------

.. list-table::
   :widths: 28 42 30
   :header-rows: 1

   * - Asset
     - Purpose
     - Download
   * - AMReX Framework Core Prompt
     - Layer 1 prompt for AMReX framework concepts and evidence boundaries.
     - :download:`Prompt <_static/AMReX_Expert_GemReport_Layer1_prompt.md>`
   * - AMReX Framework Core Report
     - Layer 1 compiled report used by the assistant.
     - :download:`Report <_static/AMReX_Expert_GemReport_Layer1-2.md>`
   * - Portable HPC / Numerics Prompt
     - Layer 2 prompt for reproducibility, numerics, and debugging policy.
     - :download:`Prompt <_static/AMReX_Expert_GemReport_Layer2_prompt.md>`
   * - Portable HPC / Numerics Report
     - Layer 2 compiled report used by the assistant.
     - :download:`Report <_static/AMReX_Expert_GemReport_Layer1-2.md>`
   * - Layer 3 Prompt A
     - Generic Layer 3 prompt optimized for ChatGPT.
     - :download:`Prompt <_static/Generic_Layer3_prompt_A.md>`
   * - Layer 3 Prompt B
     - Generic Layer 3 prompt optimized for Gemini.
     - :download:`Prompt <_static/Generic_Layer3_prompt_B.md>`
   * - Layer 3 Prompt C
     - Generic Layer 3 prompt variant from the tutorial set.
     - :download:`Prompt <_static/Generic_Layer3_prompt_C.md>`
   * - Layer 3 Prompt D
     - Generic Layer 3 prompt variant from the tutorial set.
     - :download:`Prompt <_static/Generic_Layer3_prompt_D.md>`
   * - ERF Layer 3 Generation Prompt
     - ERF-specific Layer 3 generation prompt.
     - :download:`Prompt <_static/ERF_Expert_GemReport_Layer3_prompt.md>`
   * - ERF Layer 3 Context Report
     - ERF-specific Layer 3 context report.
     - :download:`Report <_static/ERF_Expert_GemReport_Layer3.md>`
   * - ERF Layer 3 Instruction Prompt
     - ERF-specific Layer 3 instruction prompt.
     - :download:`Prompt <_static/ERF_Expert_GemReport_Layer3_instruction.md>`

Suggested setup
---------------

1. Open the ERF Specialized Assistant if you want a ready-made starting point.
2. If you want to build your own assistant, use the ERF Layer 3 context report and instruction prompt.
3. If you need to extend the assistant for a new ERF-specific workflow, update the prompt and context files together.
4. If you are writing a new simulation setup, use the assistant in inputs-authoring mode and share the scenario, domain, and hardware details first.

This keeps the assistant anchored to the same ERF sources used to generate the existing profile.
