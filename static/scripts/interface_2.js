// generate interface 
function generateInterface(listFile) {

    $("#interfaceBtn").css("display", "none")

    s = `
    <form id="box" action="http://proteinformatics.uni-leipzig.de/mutation_explorer/submit" method="post" enctype="multipart/form-data">

        Select alignment <br>
        
        <select id='fileDropdownConv' name='alignment_link' onchange='changeSeqSelectionConv()'></select>
        <br>
        <br>
        <br>

        Upload a PDB <br>

        <input id="fileConvText" name="file_conv_text" type="text"> <br>
        <input id="fileConv" name="file_conv" type="file" accept=".pdb">
        <br>
        <br>


        Upload a second PDB (optional, structures get aligned/superimposed) <br>

        <input id="fileSuperText" name="file_super_text" type="text"> <br>
        <input id="fileSuper" name="file_super" type="file" accept=".pdb">
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
                fsplit = `${files[i]}`.split("/");
                filename = fsplit[fsplit.length - 1];
                $("#fileDropdownConv").append(`<option value='${files[i]}'>${filename}</option>`);
            }
        });



    //changeSeqSelectionConv()

}








/*
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
*/
