function generateInterface(listFile) {

    $("#interfaceBtn").css("display", "none")

    s = `
    <form id="box" action="https://proteinformatics.uni-leipzig.de/voronoia/test" method="post">
        Alignment <br>
        
        <select id='fileDropdown' onchange='changeSeqSelection()'></select>
        <br>
        <br>
        <br>

        Upload 1st PDB <br>

        <select id="seqA"></select> <br>
        <br>

        <input id="pdbA" name="pdb_a" type="text" placeholder="pdb id" onchange='clearInput("#fileA")'> <br>
        <input id="fileA" name="file_a" type="file" accept=".pdb" onchange='clearInput("#pdbA")'>
        <br>
        <br>
        <br>

        Upload 2nd PDB (optional) <br>

        <select id="seqB"></select> <br>
        <br>

        <input id="pdbB" name="pdb_b" type="text" placeholder="pdb id", onchange='checkSuperimpose("#fileB")'> <br>
        <input id="fileB" name="file_b" type="file" accept=".pdb" onchange='checkSuperimpose("#pdbB")'>
        <br>
        <br>

        <input id="superimpose" type="checkbox" disabled style="margin-left: 2vw"> <span id="superimposeLabel" style="opacity: 0.5">superimpose PDBs</span>
        <br>
        <br>
        <br>

        <input id="coloring" type="checkbox"> color by conservation
        <br>
        <br>
        <br>
        <br>

        <input type="submit">
    </form>
    `;
    $("#box").append(s);


    // add alignments to dropdown
    fetch(listFile)
        .then(response => response.text())
        .then(text => text.split("\n"))
        .then(files => {
            for (let i = 0; i < files.length; i++) {
                $("#fileDropdown").append(`<option value='${files[i]}'>${files[i]}</option>`);
            }
        });



    changeSeqSelection()

}

function clearInput(el) {
    $(el).val("");
} 

function checkSuperimpose(el) {
    clearInput(el);
    if($("#pdbB").val() == "" && $("#fileB").val() == "") {
        $("#superimpose").prop("checked", false);
        $("#superimpose").attr("disabled", true);
        $("#superimposeLabel").css("opacity", "0.5");
    } else {
        $("#superimpose").attr("disabled", false);
        $("#superimposeLabel").css("opacity", "1");
    }
}


function changeSeqSelection() {
    file = $("#fileDropdown").children("option:selected").val();
    fetch(file)
        .then(response => response.text())
        //.then(text => console.log(text))
        .then(text => text.split("\n"))
        .then(lines => [lines[2].split(/[ ]+/)[0], lines[3].split(/[ ]+/)[0]])
        .then(arr => {
            $("#seqA").children().remove();
            $("#seqB").children().remove();
            for (let i = 0; i < arr.length; i++) {
                $("#seqA").append(`<option value='${arr[i]}'>${arr[i]}</option>`);
                $("#seqB").append(`<option value='${arr[i]}'>${arr[i]}</option>`);
            }
        });
}
