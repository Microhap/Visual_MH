/*
VisualMH.js - Web Program for Microhalotype Analysis with STRait Razor
Dept. of Forensic Medicine, Yonsei University College of Medicine
2021.3.28 ~ 31. 4.6, 5.5, 2022.2.14, 24 ~ 25, 7.5 ~ 18. 2023. 9. 1. 2025. 1. 14.
*/

"use strict";

var refInfo = {};

// http://www.gisdeveloper.co.kr/?p=5564
function getRefInfo(fileBlob) {
    const reader = new FileReader();

    refInfo = {};
    reader.onload = function () {
        let markerno = 0;
        let errInfo = "";

        const refLines = reader.result.split('\n');
        for (const refLine of refLines) {
            if (refLine.length < 25 || refLine.toLowerCase().indexOf("chr") < 0) continue;

            markerno += 1;
            let refField = refLine.trim().split(/\s+/);
            let seqstart = Number(refField[refField.length - 2]);
            let seqlen   = refField[refField.length - 1].length;

            let index = 0;
            let MHinfo = [];
            while (refField[index + 2].indexOf(':') > 0) {
                let MHfield = refField[index + 2].split(':');
                let snppos = Number(MHfield[1]);
                if (snppos >= seqstart && snppos < seqstart + seqlen) {
                    MHinfo.push([MHfield[0], snppos]);
                } else {
                    errInfo = errInfo + `Warning: ${refField[index + 2]} is out of ${refField[0]} range!`;
                    errInfo += '\n';
                }
                index += 1;
            }
            refInfo[refField[0]] = [refField[1], MHinfo, seqstart, refField[refField.length - 1]];
        }
        if (errInfo) alert(errInfo);

        if (markerno > 0) {
            markermsg.innerText = `${markerno} markers: Loading completed!`; // not work on IE
        } else {
            markermsg.innerText = "No available marker information!";
        }
        if (seqLines.length) getReadInfo();
    }

    if (fileBlob) {
        reader.readAsText(fileBlob, /* optional */ "euc-kr");
    } else {
        markermsg.innerText = "No information was loaded...";
        if (seqLines.length) getReadInfo();
    }
}

var seqLines = [];
var readInfo = [];
var infoLines = [];

function getSeqInfo(fileBlob) {
    const reader = new FileReader();

    seqLines = [];
    reader.onload = function () {
        seqLines = reader.result.split('\n');
        getReadInfo();
    }

    if (fileBlob) {
        reader.readAsText(fileBlob, /* optional */ "euc-kr");
    } else {
        resultmsg.innerText = "No sequence was analyzed...";
        infoLines = [];
        // resultlist.innerText = infoLines;
        showResult();
    }
}

var minreadno, noisecut, homopolyerr, allelecover, minormarked;

function getOption() {
    minreadno = Number(minreadnoval.value);
    noisecut  = Number(noisecutval.value);
    homopolyerr = Number(homopolyval.value);
    allelecover = Number(alleleCRval.value);
    minormarked = alleleminor.checked;

    if (seqLines.length) getReadInfo();
}

function getReadInfo() {
    readInfo = [];
    let currmarker = "";
    let markerno = 0;
    let coverage = 0;
    infoLines = [];

    for (let i in seqLines) {
        if (seqLines[i].length < 30 || seqLines[i].indexOf("bases") < 0) continue;

        let seqFields = seqLines[i].trim().split(/\s+/);
        let markerallele = seqFields[0].split(':');
        if (currmarker && currmarker !== markerallele[0]) {
            markerno += 1;
            parseSeq(currmarker, coverage);
            readInfo = [];
            coverage = 0;
        }

        let readcount = Number(seqFields[4]) + Number(seqFields[5]);
        if (readcount >= minreadno && seqFields[1] !== "0") // Check Below Theshold
            readInfo.push([markerallele[1], seqFields[3], readcount]);
        coverage = coverage + readcount;
        currmarker = markerallele[0];
    }

    if (readInfo.length > 0) {
        markerno += 1;
        parseSeq(currmarker, coverage);
    }

    if (markerno > 0) {
        resultmsg.innerText = `${markerno} markers: Analysis completed!`; // not work on IE
    } else {
        resultmsg.innerText = "No valid marker!";
    }
    // resultlist.innerText = infoLines;
    showResult();
}

const rsNumberOn = true;

function parseSeq(marker, coverage) {
    let snpInfo;

    if (marker in refInfo) {
        let refvalue = refInfo[marker];
        let refseq = refvalue[3].toUpperCase();
        let reflen = refseq.length;
        let SNPpos = refvalue[1].map(MHinfo => MHinfo[1] - refvalue[2]);
        // for (const MHinfo of refvalue[1]) SNPpos.push(MHinfo[1] - refvalue[2]);

        let validread = 0;
        for (const readLine of readInfo)        
            if (readLine[2] / coverage > noisecut) validread = validread + readLine[2];

        let alleleno = 0;
        if (validread > 0) {
            for (const [_, sequence, readcount] of readInfo) {
                let validratio = readcount / validread;
                let grossratio = readcount / coverage;

                if (grossratio <= noisecut) continue;  // Check Noise Read
                let hpolyflag = isHomopolymer(sequence);
                if (hpolyflag && validratio <= homopolyerr) continue;
                if (validratio <= allelecover && !minormarked) continue;

                alleleno += 1;
                let alignstr = alignseq(refseq, sequence);

                let rsinfo = "";
                let refsnp = "";
                let seqsnp = "";
                snpInfo = `${refvalue[0]}:${refvalue[1][0][1]}`;
                SNPpos.forEach(function (pos, i) {
                    rsinfo = rsinfo + ' ' + `${refvalue[1][i][0]}:`;
                    if (refseq[pos] === alignstr[pos]) {
                        rsinfo = rsinfo + refseq[pos] + "  ";
                    } else {
                        rsinfo = rsinfo + refseq[pos] + "> ";
                    }
                    refsnp = refsnp + refseq[pos];  // <-> refvalue[3][pos]
                    seqsnp = seqsnp + alignstr[pos];
                })
                if (rsNumberOn) {
                    snpInfo = snpInfo + '\t' + rsinfo.replace(/^\s+/, '') + '\t' + seqsnp;
                } else {
                    snpInfo = snpInfo + '\t' + refsnp + '\t' + seqsnp;
                }

                let extInfo = "";
                if (indelpos >= 0) {
                    let indelflag = (sequence.length < reflen) ? '-' : '+';
					if (indelflag === '+' || !SNPpos.includes(indelpos)) {
						extInfo = extInfo + ' ' + `${refvalue[2] + indelpos}${indelflag}[${indelstr}]`;
						if (hpolyflag) extInfo += '@';
					}
                }
                for (let i = 0; i < reflen; i++) {
                    if (diffpos.includes(i) && !SNPpos.includes(i))  // <-> i in diffpos and i not in SNPpos
                        extInfo = extInfo + ' ' + `${refvalue[2] + i}` + refseq[i] + '>' + alignstr[i];
                }
                snpInfo = snpInfo + '\t' + extInfo.replace(/^\s+/, '');  // ltrim()

                let minorflag = (validratio <= allelecover) ? '<' : ' ';  // Check Read Balance
                snpInfo = `${marker}:${alleleno}\t${snpInfo}\t${readcount}\t${validratio.toFixed(3)}\t${grossratio.toFixed(3)}\t${minorflag}`;
                infoLines.push(snpInfo);
            }
        } else {
            snpInfo = `${marker}:0.0\t${coverage}\tNo valid read!`;
            infoLines.push(snpInfo);
        }
    } else {
        snpInfo = `${marker}\t${coverage}\tNo Ref. Info.`;
        infoLines.push(snpInfo);
    }
    infoLines.push("");
}

function isHomopolymer(seqstr) {
    const polyNo = 8;

    return seqstr.includes('A'.repeat(polyNo)) || seqstr.includes('T'.repeat(polyNo)) || 
           seqstr.includes('G'.repeat(polyNo)) || seqstr.includes('C'.repeat(polyNo)) ||
           seqstr.includes('AT'.repeat(5));
}

var indelpos, indelstr;
var diffpos = [];

function alignseq(refstr, seqstr) {
    indelpos = -1;
    indelstr = "";
    let gabchr = '-';

    let reflen = refstr.length;
    let seqlen = seqstr.length;
    let indellen = Math.abs(reflen - seqlen);
    const score = new Array(seqlen + 1).fill(0);

    let alignstr = seqstr;
    if (reflen > seqlen) {
        for (let i = 0; i < seqlen + 1; i++) {
            alignstr = seqstr.slice(0, i) + gabchr.repeat(indellen) + seqstr.slice(i);
            for (let j = 0; j < reflen; j++)
                if (refstr[j] === alignstr[j]) score[i] += 1;
                else if (alignstr[j] !== gabchr) score[i] -= 1;
        }
        indelpos = score.indexOf(Math.max(...score));
        indelstr = refstr.slice(indelpos, indelpos + indellen);
        alignstr = seqstr.slice(0, indelpos) + gabchr.repeat(indellen) + seqstr.slice(indelpos);
    } else
    if (reflen < seqlen) {
        for (let i = 0; i < seqlen - indellen + 1; i++) {
            alignstr = seqstr.slice(0, i) + seqstr.slice(i + indellen);
            for (let j = 0; j < reflen; j++)
                if (refstr[j] === alignstr[j]) score[i] += 1;
                else score[i] -= 1;
        }
        indelpos = score.indexOf(Math.max(...score));
        indelstr = seqstr.slice(indelpos, indelpos + indellen);
        alignstr = seqstr.slice(0, indelpos) + seqstr.slice(indelpos + indellen);
    }

    diffpos = [];
    for (let i = 0; i < reflen; i++)
        if (refstr[i] !== alignstr[i] && alignstr[i] !== gabchr)
            diffpos.push(i);

    return alignstr;
}

function showResult() {
    const tbodyList = document.querySelector("tbody");
    while (tbodyList.hasChildNodes())
        tbodyList.firstChild.remove();  // tbodyList.removeChild(tbodyList.firstChild);

    const template = document.querySelector("template");
    for (let i = 0; i < infoLines.length - 1; i++) {
        const tempNode = document.importNode(template.content, true);
        const trRow  = tempNode.querySelector("tr");
        const tdCols = tempNode.querySelectorAll("td");

        if (infoLines[i]) {
            trRow.style.background = "#FAFAFA";
            // if (infoLines[i].includes('@')) trRow.style.background = "mintcream";
            if (infoLines[i].includes('<')) trRow.style.background = "lightyellow";

            let lineField = infoLines[i].split('\t');
            for (let j in lineField) {
                if (lineField[j].includes("No")) tdCols[j].colSpan = 3;
                tdCols[j].textContent = lineField[j];
            }
        }
        tbodyList.append(tempNode);  // tbodyList.appendChild(tempNode);
    }
}


// http://www.gisdeveloper.co.kr/?p=5564
function saveToFile_Chrome(fileName, content) {
    const blob = new Blob([content], { type: "text/plain", endings: "native" });
    const objURL = window.URL.createObjectURL(blob);

    // 이전에 생성된 메모리 해제
    if (window.__Xr_objURL_forCreatingFile__)
        window.URL.revokeObjectURL(window.__Xr_objURL_forCreatingFile__);
    window.__Xr_objURL_forCreatingFile__ = objURL;

    const link = document.createElement("a");
    link.download = fileName;
    link.href = objURL;
    link.click();
}


// window.onload = function() {
window.addEventListener("load", function () {
    title.onclick = () => location.reload();

    const resetOption = document.getElementById("btn-reset");
    const optionValue = document.getElementById("option-info");
    const numberBoxes = document.querySelectorAll(".number-only");

    resetOption.onclick = function () {
        minreadnoval.value = "100";
        noisecutval.value = "0.05";
        homopolyval.value = "0.15";  // confirmed on 3/30!
        alleleCRval.value = "0.25";  // 1/3 of major read
        alleleminor.checked = true;

        getOption();
    };

    optionValue.onchange = function (e) {
        if (e.target.id === "minreadnoval") {
            if (Number(e.target.value) < Number(minreadnoval.min)) {
                e.target.value = minreadnoval.min;  // 10
                minreadnoval.focus();
            }
        } else {
            e.target.value = e.target.value.replace(/[^0-9\.]/g, "");
        }

        getOption();
    }

    numberBoxes.forEach (function (numberBox) {
        // https://www.3rabbitz.com/blog_ko/a6068a3f953de9d8
        numberBox.onkeydown = function (e) {
            const acceptCode = [9, 36, 35, 37, 39, 8, 46];

            if ((e.key >= '0' && e.key <= '9') || e.key === '.') return;
            if (acceptCode.includes(e.keyCode) || e.ctrlKey || e.altKey) return;
            if (e.keyCode >= 112 && e.keyCode <= 123 || e.keyCode === 13) return;
            e.preventDefault();
        }
    })

    const currRefInfo = document.getElementById("curr-marker");
    const currSeqFile = document.getElementById("curr-sequence");
    const btnDownload = document.getElementById("btn-download");

    currRefInfo.accept = "text/plain"; // 확장자가 xxx, yyy 일때, ".xxx, .yyy"
    currSeqFile.accept = "text/plain"; // 확장자가 xxx, yyy 일때, ".xxx, .yyy"

    currRefInfo.onchange = (e) => getRefInfo(e.target.files[0]);
    currSeqFile.onchange = (e) => getSeqInfo(e.target.files[0]);

    btnDownload.onclick = function () {
        if (!currSeqFile.files[0]) return;

        let fileName = currSeqFile.files[0].name;
        if (fileName.indexOf("_Seq.") > 0) {
            fileName = fileName.replace("_Seq.", "_MH.");
        } else {
            fileName = fileName.replace(".txt", "_MH.txt");
        }
        if (infoLines.length) saveToFile_Chrome(fileName, infoLines.join('\n'));
    }

    resetOption.click();
})
