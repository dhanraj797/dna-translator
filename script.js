(function () {
  "use strict";
  // Genetic code, codons ordered U, C, A, G
  var AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  var B = "UCAG", CODON = {};
  for (var i = 0; i < 64; i++) CODON[B[i >> 4] + B[(i >> 2) & 3] + B[i & 3]] = AA[i];
  var NAMES = {A:"Ala",R:"Arg",N:"Asn",D:"Asp",C:"Cys",E:"Glu",Q:"Gln",G:"Gly",H:"His",I:"Ile",L:"Leu",K:"Lys",M:"Met",F:"Phe",P:"Pro",S:"Ser",T:"Thr",W:"Trp",Y:"Tyr",V:"Val","*":"Stop"};
  var CLS = {"*": "stop"};
  "GAVLIMFWP".split("").forEach(function (c) { CLS[c] = "non"; });
  "STCYNQ".split("").forEach(function (c) { CLS[c] = "pol"; });
  "KRH".split("").forEach(function (c) { CLS[c] = "pos"; });
  "DE".split("").forEach(function (c) { CLS[c] = "neg"; });

  var SAMPLE = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG";
  var $ = function (id) { return document.getElementById(id); };
  var current = { mrna: "", protein: "" };

  function esc(s) { return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;"); }

  // Step 1: validation. Drops FASTA headers, whitespace and digits, then flags anything not A, C, G or T.
  function validate(raw) {
    var text = raw.split(/\r?\n/).filter(function (l) { return l.trim().charAt(0) !== ">"; }).join("");
    var seq = text.toUpperCase().replace(/[\s\d]/g, "");
    var invalid = [];
    for (var i = 0; i < seq.length; i++) if ("ACGT".indexOf(seq[i]) < 0) invalid.push({ ch: seq[i], pos: i + 1 });
    return { seq: seq, invalid: invalid };
  }
  function revcompDNA(s) {
    var m = { A: "T", T: "A", C: "G", G: "C" }, o = "";
    for (var i = s.length - 1; i >= 0; i--) o += m[s[i]];
    return o;
  }
  function revcompRNA(s) {
    var m = { A: "U", U: "A", C: "G", G: "C" }, o = "";
    for (var i = s.length - 1; i >= 0; i--) o += m[s[i]];
    return o;
  }
  // Step 2: transcription. mRNA matches the coding strand with T replaced by U.
  function transcribe(dna, strand) {
    var coding = strand === "template" ? revcompDNA(dna) : dna;
    return coding.replace(/T/g, "U");
  }
  // Step 3: translation.
  function codonsOf(mrna, frame) {
    var s = frame >= 3 ? revcompRNA(mrna) : mrna, out = [];
    for (var i = frame % 3; i + 3 <= s.length; i += 3) { var c = s.substr(i, 3); out.push([c, CODON[c]]); }
    return out;
  }
  function trim(list, atg, stp) {
    var start = 0, out;
    if (atg) {
      start = -1;
      for (var i = 0; i < list.length; i++) if (list[i][0] === "AUG") { start = i; break; }
      if (start < 0) return [];
    }
    out = list.slice(start);
    if (stp) for (var j = 0; j < out.length; j++) if (out[j][1] === "*") { out = out.slice(0, j + 1); break; }
    return out;
  }
  function clearOutputs(msg) {
    current = { mrna: "", protein: "" };
    $("mrnaOut").textContent = ""; $("protein").textContent = msg || "";
    $("chips").innerHTML = ""; $("frames").innerHTML = ""; $("stats").innerHTML = "";
  }
  function setStatus(msg, bad) { var s = $("status"); s.textContent = msg; s.className = "status" + (bad ? " bad" : ""); }

  function render() {
    var raw = $("seq").value;
    if (!raw.trim()) { $("dnaOut").textContent = ""; setStatus("Paste a DNA sequence to begin."); clearOutputs("Your protein will appear here."); return; }
    var v = validate(raw);

    if (v.invalid.length) {
      var html = "";
      for (var k = 0; k < v.seq.length; k++) html += "ACGT".indexOf(v.seq[k]) < 0 ? "<mark>" + esc(v.seq[k]) + "</mark>" : v.seq[k];
      $("dnaOut").innerHTML = html;
      var hasU = v.invalid.some(function (x) { return x.ch === "U"; });
      setStatus(v.invalid.length + " invalid character" + (v.invalid.length > 1 ? "s" : "") + ", first at position " + v.invalid[0].pos +
        ". DNA accepts only A, C, G and T." + (hasU ? " U belongs in RNA; enter the DNA version with T." : ""), true);
      clearOutputs("Fix the highlighted characters to continue.");
      return;
    }

    var dna = v.seq, frame = +$("frame").value;
    $("dnaOut").textContent = dna;
    setStatus("Valid DNA, " + dna.length + " bases." + (dna.length % 3 ? " Length is not a multiple of 3; trailing bases are ignored." : ""));

    var mrna = transcribe(dna, $("strand").value);
    var list = trim(codonsOf(mrna, frame), $("atg").checked, $("stp").checked);
    var prot = list.map(function (x) { return x[1]; }).join("");
    current = { mrna: mrna, protein: prot.replace(/\*/g, "") };

    $("mrnaOut").textContent = mrna;
    $("protein").textContent = prot || "No protein found. Try another frame, or turn off Start at first AUG.";
    var gc = Math.round(100 * dna.replace(/[AT]/g, "").length / dna.length);
    $("stats").innerHTML = "<span><b>" + dna.length + "</b>bases</span><span><b>" + gc + "%</b>GC</span><span><b>" + current.protein.length + "</b>amino acids</span>";

    var h = "";
    list.forEach(function (x, i) {
      h += '<div class="chip ' + CLS[x[1]] + (i === 0 && x[0] === "AUG" && $("atg").checked ? " start" : "") + '"><div class="cd">' + x[0] +
        '</div><div class="aa">' + x[1] + '</div><div class="nm">' + NAMES[x[1]] + "</div></div>";
    });
    $("chips").innerHTML = h;

    var fh = "";
    for (var f = 0; f < 6; f++) {
      var p = codonsOf(mrna, f).map(function (x) { return x[1]; }).join("");
      fh += '<button type="button" data-f="' + f + '" class="' + (f === frame ? "on" : "") + '"><b>' + (f < 3 ? "+" : "\u2212") + (f % 3 + 1) + "</b><span>" + (p || "-") + "</span></button>";
    }
    $("frames").innerHTML = fh;
  }

  $("seq").addEventListener("input", render);
  ["strand", "frame", "atg", "stp"].forEach(function (id) { $(id).addEventListener("change", render); });
  $("sample").addEventListener("click", function () { $("seq").value = SAMPLE; $("strand").value = "coding"; $("frame").value = "0"; render(); });
  $("clear").addEventListener("click", function () { $("seq").value = ""; render(); $("seq").focus(); });
  $("frames").addEventListener("click", function (e) {
    var b = e.target.closest("button"); if (!b) return; $("frame").value = b.getAttribute("data-f"); render();
  });
  document.addEventListener("click", function (e) {
    var b = e.target.closest("[data-copy]"); if (!b) return;
    var label = b.textContent, text = current[b.getAttribute("data-copy")];
    function done(m) { b.textContent = m; setTimeout(function () { b.textContent = label; }, 1500); }
    if (navigator.clipboard) navigator.clipboard.writeText(text).then(function () { done("Copied"); }, function () { done("Copy failed"); });
    else done("Copy failed");
  });
  render();
})();
