import{_ as a,o as e,c as n,aA as i}from"./chunks/framework.DuD_VnaT.js";const k=JSON.parse('{"title":"Realignment","description":"","frontmatter":{},"headers":[],"relativePath":"demos/erp/realign.md","filePath":"demos/erp/realign.md","lastUpdated":null}'),t={name:"demos/erp/realign.md"};function l(p,s,o,r,h,c){return e(),n("div",null,[...s[0]||(s[0]=[i(`<h1 id="Realignment" tabindex="-1">Realignment <a class="header-anchor" href="#Realignment" aria-label="Permalink to &quot;Realignment {#Realignment}&quot;">​</a></h1><p>This demo shows how to realign stimulus-locked epochs to a different time point, such as response time, for response-locked ERP analysis.</p><h3 id="When-to-Use-Realignment" tabindex="-1">When to Use Realignment <a class="header-anchor" href="#When-to-Use-Realignment" aria-label="Permalink to &quot;When to Use Realignment {#When-to-Use-Realignment}&quot;">​</a></h3><p>Stimulus-locked epochs are aligned to stimulus onset (t=0). Realignment shifts t=0 to a different event — typically the participant&#39;s response — so that you can study activity time-locked to that event instead.</p><p>Common use cases:</p><ul><li><p><strong>Response-locked ERPs</strong> — study motor preparation relative to button press</p></li><li><p><strong>Saccade-locked ERPs</strong> — study activity relative to eye movement onset</p></li><li><p><strong>Any event-locked analysis</strong> — realign to any column in the epoch data</p></li></ul><h3 id="How-it-Works" tabindex="-1">How it Works <a class="header-anchor" href="#How-it-Works" aria-label="Permalink to &quot;How it Works {#How-it-Works}&quot;">​</a></h3><ol><li><p>Each epoch&#39;s time vector is shifted so that the realignment value becomes t=0</p></li><li><p>All epochs are cropped to the common time interval that is valid across all trials</p></li><li><p>A uniform time vector is regenerated to ensure consistency</p></li></ol><h3 id="Key-Functions" tabindex="-1">Key Functions <a class="header-anchor" href="#Key-Functions" aria-label="Permalink to &quot;Key Functions {#Key-Functions}&quot;">​</a></h3><table tabindex="0"><thead><tr><th style="text-align:right;">Function</th><th style="text-align:right;">Purpose</th></tr></thead><tbody><tr><td style="text-align:right;"><code>realign!(epochs, :rt)</code></td><td style="text-align:right;">Realign in place (mutating)</td></tr><tr><td style="text-align:right;"><code>realign(epochs, :rt)</code></td><td style="text-align:right;">Return a realigned copy</td></tr><tr><td style="text-align:right;"><code>realign(file_pattern, :rt)</code></td><td style="text-align:right;">Batch realign across participants</td></tr></tbody></table><h2 id="Workflow-Summary" tabindex="-1">Workflow Summary <a class="header-anchor" href="#Workflow-Summary" aria-label="Permalink to &quot;Workflow Summary {#Workflow-Summary}&quot;">​</a></h2><h3 id="Single-Participant-Realignment" tabindex="-1">Single-Participant Realignment <a class="header-anchor" href="#Single-Participant-Realignment" aria-label="Permalink to &quot;Single-Participant Realignment {#Single-Participant-Realignment}&quot;">​</a></h3><ul><li>Realign epochs to response time column</li></ul><h3 id="Batch-Realignment" tabindex="-1">Batch Realignment <a class="header-anchor" href="#Batch-Realignment" aria-label="Permalink to &quot;Batch Realignment {#Batch-Realignment}&quot;">​</a></h3><ul><li>Process all participant files in a directory</li></ul><h3 id="Typical-Pipeline" tabindex="-1">Typical Pipeline <a class="header-anchor" href="#Typical-Pipeline" aria-label="Permalink to &quot;Typical Pipeline {#Typical-Pipeline}&quot;">​</a></h3><ul><li>Extract stimulus-locked epochs → realign to RT → average → LRP → jackknife</li></ul><h2 id="Code-Examples" tabindex="-1">Code Examples <a class="header-anchor" href="#Code-Examples" aria-label="Permalink to &quot;Code Examples {#Code-Examples}&quot;">​</a></h2><details class="details custom-block"><summary>Show Code</summary><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Demo: Response-Locked Realignment</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Shows how to realign stimulus-locked epochs to a different time point</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># (e.g., response time) for response-locked ERP analysis.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> EegFun</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#######################################################################</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># SINGLE-PARTICIPANT REALIGNMENT</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#######################################################################</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Load stimulus-locked epoched data</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># epochs = EegFun.read_data(&quot;participant_1_epochs.jld2&quot;)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Realign to response times stored in the :rt column</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># EegFun.realign!(epochs, :rt)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Non-mutating version (returns a new copy)</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># realigned = EegFun.realign(epochs, :rt)</span></span>
<span class="line"></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#######################################################################</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># BATCH REALIGNMENT</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#######################################################################</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Process all participant epoch files in a directory.</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Each epoch must have a column with the realignment time</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># (e.g., response time relative to stimulus onset).</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># EegFun.realign(&quot;epochs_cleaned&quot;, :rt)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Specify input directory</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># EegFun.realign(&quot;epochs_cleaned&quot;, :rt, input_dir = &quot;/path/to/epochs/&quot;)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Specific participants</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># EegFun.realign(&quot;epochs_cleaned&quot;, :rt, participant_selection = EegFun.participants([1, 2, 3]))</span></span>
<span class="line"></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#######################################################################</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># TYPICAL WORKFLOW</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#######################################################################</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Response-locked LRP analysis:</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Extract stimulus-locked epochs with RT column</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Realign to response time:</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#      realign(&quot;epochs_cleaned&quot;, :rt)</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Average to ERPs:</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#      average_epochs in realigned directory</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Calculate LRP:</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#      lrp(&quot;realigned_erps&quot;, [(1, 2)])</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Jackknife average:</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#      jackknife_average(&quot;lrp&quot;)</span></span></code></pre></div></details><h2 id="See-Also" tabindex="-1">See Also <a class="header-anchor" href="#See-Also" aria-label="Permalink to &quot;See Also {#See-Also}&quot;">​</a></h2><ul><li><a href="./../../reference/">API Reference</a></li></ul>`,21)])])}const g=a(t,[["render",l]]);export{k as __pageData,g as default};
