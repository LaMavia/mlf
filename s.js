const http = require("https"); // or 'https' for https:// URLs
const fs = require("fs");

const file = fs.createWriteStream("file");
const request = http.get(
  "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28trmd%29",
  function (response) {
    response.pipe(file);

    // after download completed close filestream
    file.on("finish", () => {
      file.close();
      console.log("Download Completed");
    });
  },
);
