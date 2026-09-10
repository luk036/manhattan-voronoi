import {generateL1Voronoi, generateL2Voronoi} from "./src/voronoi.js";

const width = 400;
const height = 400;
const svgNs = "http://www.w3.org/2000/svg";

function randomNormal(sharpness){
    return new Array(sharpness).fill(0).map(() => Math.random()).reduce((c, e) => c + e, 0) / sharpness;
}

const params = (new URL(document.location)).searchParams;
const numPoints = parseInt(params.get("points")) || 64;

const seen = new Set();
const sites = new Array(numPoints).fill(0)
    .map(() => [Math.floor(randomNormal(2) * width), Math.floor(randomNormal(2) * height)])
    .filter(point => {
        const key = point.join(",");
        if(seen.has(key)){
            return false;
        }
        seen.add(key);
        return true;
    });

document.getElementById("points").textContent = JSON.stringify(
    sites.slice().sort((a, b) => a[0] !== b[0] ? a[0] - b[0] : a[1] - b[1]), null);

["l1-diagram", "l1-delaunay", "l2-diagram", "l2-delaunay"].forEach(id => {
    const svg = document.getElementById(id);
    svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
});

function drawCells(svg, cells, fill){
    cells.forEach(cell => {
        const path = document.createElementNS(svgNs, "path");
        path.setAttribute("d", cell.d);
        path.setAttribute("class", "polygon");
        path.style.fill = fill;
        path.style.stroke = "#000";
        path.style.strokeWidth = "1px";
        svg.appendChild(path);
    });
}

function drawSites(svg, points){
    points.forEach(point => {
        const circle = document.createElementNS(svgNs, "circle");
        circle.setAttribute("cx", point[0]);
        circle.setAttribute("cy", point[1]);
        circle.setAttribute("r", 1.5);
        circle.style.fill = "#000";
        svg.appendChild(circle);
    });
}

function drawDelaunay(svg, cells, color){
    const drawn = new Set();
    cells.forEach(cell => {
        cell.neighbors.forEach(neighbor => {
            const key = [cell.site.join(","), neighbor.join(",")].sort().join("|");
            if(drawn.has(key)){
                return;
            }
            drawn.add(key);
            const line = document.createElementNS(svgNs, "line");
            line.setAttribute("x1", cell.site[0]);
            line.setAttribute("y1", cell.site[1]);
            line.setAttribute("x2", neighbor[0]);
            line.setAttribute("y2", neighbor[1]);
            line.style.stroke = color;
            line.style.strokeWidth = "1px";
            svg.appendChild(line);
        });
    });
    drawSites(svg, cells.map(cell => cell.site));
}

const l1Cells = generateL1Voronoi(sites, width, height, true);
drawCells(document.getElementById("l1-diagram"), l1Cells, "#e3f2fd");
drawSites(document.getElementById("l1-diagram"), l1Cells.map(cell => cell.site));
drawDelaunay(document.getElementById("l1-delaunay"), l1Cells, "#1565c0");

const l2Cells = generateL2Voronoi(sites, width, height, true);
drawCells(document.getElementById("l2-diagram"), l2Cells, "#fff3e0");
drawSites(document.getElementById("l2-diagram"), l2Cells.map(cell => cell.site));
drawDelaunay(document.getElementById("l2-delaunay"), l2Cells, "#c2185b");
