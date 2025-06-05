// This script automatically generates a table of contents dynamically. 
// It is a slightly modified version of a script from: 
// http://blog.magnetiq.com/post/497600148/automatic-table-of-contents-generation
// Currently this will not grab a header to put in the table of contents if
// it has any embeded html tags in the header text. It's just a case of 
// fixing the regular experession to get this to work properly
// Got the negating regular expression to allow other tags from:
// http://stackoverflow.com/questions/977251/regular-expressions-and-negating-a-whole-character-group
// Gives (?:(?!<h).*)+

window.onload = function () {
        var toc = "";
        var level = 0;

        <!-- original regular expression: <h([\d]) .*>([^<]+)<\/h([\d])> -->
        document.getElementById("contents").innerHTML =
                document.getElementById("contents").innerHTML.replace(
                      /<h([\d]) .*>([^<]+)<\/h([\d])>/gi, 
                        function (str, openLevel, titleText, closeLevel) {
                                if (openLevel != closeLevel) {
                                        return str;
                                }

                                if (openLevel > level) {
                                toc += (new Array(openLevel - level + 1)).join("<ul>");
                                } else if (openLevel < level) {
                                        toc += (new Array(level - openLevel + 1))
                        .join("</ul>");
                                }

                                level = parseInt(openLevel);

                                var anchor = titleText.replace(/ /g, "_");
                                toc += "<li><a href=\"#" + anchor + "\">"
                    + titleText + "</a></li>";

                                return "<h" + openLevel + " class=\"numbered\" id=\"" + anchor
                    + "\">" + titleText + "</h" + closeLevel
                    + ">";
                        }
                );

        if (level) {
                toc += (new Array(level + 1)).join("</ul>");
        }

        document.getElementById("toc").innerHTML += toc;
};
