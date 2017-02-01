// Apache-Style Software License for ColorBrewer software and ColorBrewer Color
// Schemes
//
// Copyright (c) 2002 Cynthia Brewer, Mark Harrower, and The Pennsylvania State
// University.
// Copyright (c) 2015 Christoph Schulz.
//
// Licensed under the Apache License, Version 2.0 (the "License"}; you may not
// use this file except in compliance with the License. You may obtain a copy of
// the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
// WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
// License for the specific language governing permissions and limitations under
// the License.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions as source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. The end-user documentation included with the redistribution, if any, must
// include the following acknowledgment: "This product includes color
// specifications and designs developed by Cynthia Brewer
// (http://colorbrewer.org/)." Alternately, this acknowledgment may appear in the
// software itself, if and wherever such third-party acknowledgments normally
// appear.
//
// 4. The name "ColorBrewer" must not be used to endorse or promote products
// derived from this software without prior written permission. For written
// permission, please contact Cynthia Brewer at cbrewer@psu.edu.
//
// 5. Products derived from this software may not be called "ColorBrewer", nor
// may "ColorBrewer" appear in their name, without prior written permission of
// Cynthia Brewer.
#ifndef COLORBREWER_H
#define COLORBREWER_H

#include <string>
#include <initializer_list>

template <typename T>
inline std::initializer_list<T> brew(const std::string& colorName, size_t colorCount)
{
    if (colorName == "YlGn") {
        switch (colorCount) {
        case 3:
            return { "#f7fcb9", "#addd8e", "#31a354" };
        case 4:
            return { "#ffffcc", "#c2e699", "#78c679", "#238443" };
        case 5:
            return { "#ffffcc", "#c2e699", "#78c679", "#31a354", "#006837" };
        case 6:
            return { "#ffffcc", "#d9f0a3", "#addd8e", "#78c679", "#31a354", "#006837" };
        case 7:
            return { "#ffffcc", "#d9f0a3", "#addd8e", "#78c679", "#41ab5d", "#238443", "#005a32" };
        case 8:
            return { "#ffffe5", "#f7fcb9", "#d9f0a3", "#addd8e", "#78c679", "#41ab5d", "#238443", "#005a32" };
        case 9:
            return { "#ffffe5", "#f7fcb9", "#d9f0a3", "#addd8e", "#78c679", "#41ab5d", "#238443", "#006837", "#004529" };
        }
    }
    else if (colorName == "YlGnBu") {
        switch (colorCount) {
        case 3:
            return { "#edf8b1", "#7fcdbb", "#2c7fb8" };
        case 4:
            return { "#ffffcc", "#a1dab4", "#41b6c4", "#225ea8" };
        case 5:
            return { "#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494" };
        case 6:
            return { "#ffffcc", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494" };
        case 7:
            return { "#ffffcc", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84" };
        case 8:
            return { "#ffffd9", "#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84" };
        case 9:
            return { "#ffffd9", "#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58" };
        }
    }
    else if (colorName == "GnBu") {
        switch (colorCount) {
        case 3:
            return { "#e0f3db", "#a8ddb5", "#43a2ca" };
        case 4:
            return { "#f0f9e8", "#bae4bc", "#7bccc4", "#2b8cbe" };
        case 5:
            return { "#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac" };
        case 6:
            return { "#f0f9e8", "#ccebc5", "#a8ddb5", "#7bccc4", "#43a2ca", "#0868ac" };
        case 7:
            return { "#f0f9e8", "#ccebc5", "#a8ddb5", "#7bccc4", "#4eb3d3", "#2b8cbe", "#08589e" };
        case 8:
            return { "#f7fcf0", "#e0f3db", "#ccebc5", "#a8ddb5", "#7bccc4", "#4eb3d3", "#2b8cbe", "#08589e" };
        case 9:
            return { "#f7fcf0", "#e0f3db", "#ccebc5", "#a8ddb5", "#7bccc4", "#4eb3d3", "#2b8cbe", "#0868ac", "#084081" };
        }
    }
    else if (colorName == "BuGn") {
        switch (colorCount) {
        case 3:
            return { "#e5f5f9", "#99d8c9", "#2ca25f" };
        case 4:
            return { "#edf8fb", "#b2e2e2", "#66c2a4", "#238b45" };
        case 5:
            return { "#edf8fb", "#b2e2e2", "#66c2a4", "#2ca25f", "#006d2c" };
        case 6:
            return { "#edf8fb", "#ccece6", "#99d8c9", "#66c2a4", "#2ca25f", "#006d2c" };
        case 7:
            return { "#edf8fb", "#ccece6", "#99d8c9", "#66c2a4", "#41ae76", "#238b45", "#005824" };
        case 8:
            return { "#f7fcfd", "#e5f5f9", "#ccece6", "#99d8c9", "#66c2a4", "#41ae76", "#238b45", "#005824" };
        case 9:
            return { "#f7fcfd", "#e5f5f9", "#ccece6", "#99d8c9", "#66c2a4", "#41ae76", "#238b45", "#006d2c", "#00441b" };
        }
    }
    else if (colorName == "PuBuGn") {
        switch (colorCount) {
        case 3:
            return { "#ece2f0", "#a6bddb", "#1c9099" };
        case 4:
            return { "#f6eff7", "#bdc9e1", "#67a9cf", "#02818a" };
        case 5:
            return { "#f6eff7", "#bdc9e1", "#67a9cf", "#1c9099", "#016c59" };
        case 6:
            return { "#f6eff7", "#d0d1e6", "#a6bddb", "#67a9cf", "#1c9099", "#016c59" };
        case 7:
            return { "#f6eff7", "#d0d1e6", "#a6bddb", "#67a9cf", "#3690c0", "#02818a", "#016450" };
        case 8:
            return { "#fff7fb", "#ece2f0", "#d0d1e6", "#a6bddb", "#67a9cf", "#3690c0", "#02818a", "#016450" };
        case 9:
            return { "#fff7fb", "#ece2f0", "#d0d1e6", "#a6bddb", "#67a9cf", "#3690c0", "#02818a", "#016c59", "#014636" };
        }
    }
    else if (colorName == "PuBu") {
        switch (colorCount) {
        case 3:
            return { "#ece7f2", "#a6bddb", "#2b8cbe" };
        case 4:
            return { "#f1eef6", "#bdc9e1", "#74a9cf", "#0570b0" };
        case 5:
            return { "#f1eef6", "#bdc9e1", "#74a9cf", "#2b8cbe", "#045a8d" };
        case 6:
            return { "#f1eef6", "#d0d1e6", "#a6bddb", "#74a9cf", "#2b8cbe", "#045a8d" };
        case 7:
            return { "#f1eef6", "#d0d1e6", "#a6bddb", "#74a9cf", "#3690c0", "#0570b0", "#034e7b" };
        case 8:
            return { "#fff7fb", "#ece7f2", "#d0d1e6", "#a6bddb", "#74a9cf", "#3690c0", "#0570b0", "#034e7b" };
        case 9:
            return { "#fff7fb", "#ece7f2", "#d0d1e6", "#a6bddb", "#74a9cf", "#3690c0", "#0570b0", "#045a8d", "#023858" };
        }
    }
    else if (colorName == "BuPu") {
        switch (colorCount) {
        case 3:
            return { "#e0ecf4", "#9ebcda", "#8856a7" };
        case 4:
            return { "#edf8fb", "#b3cde3", "#8c96c6", "#88419d" };
        case 5:
            return { "#edf8fb", "#b3cde3", "#8c96c6", "#8856a7", "#810f7c" };
        case 6:
            return { "#edf8fb", "#bfd3e6", "#9ebcda", "#8c96c6", "#8856a7", "#810f7c" };
        case 7:
            return { "#edf8fb", "#bfd3e6", "#9ebcda", "#8c96c6", "#8c6bb1", "#88419d", "#6e016b" };
        case 8:
            return { "#f7fcfd", "#e0ecf4", "#bfd3e6", "#9ebcda", "#8c96c6", "#8c6bb1", "#88419d", "#6e016b" };
        case 9:
            return { "#f7fcfd", "#e0ecf4", "#bfd3e6", "#9ebcda", "#8c96c6", "#8c6bb1", "#88419d", "#810f7c", "#4d004b" };
        }
    }
    else if (colorName == "RdPu") {
        switch (colorCount) {
        case 3:
            return { "#fde0dd", "#fa9fb5", "#c51b8a" };
        case 4:
            return { "#feebe2", "#fbb4b9", "#f768a1", "#ae017e" };
        case 5:
            return { "#feebe2", "#fbb4b9", "#f768a1", "#c51b8a", "#7a0177" };
        case 6:
            return { "#feebe2", "#fcc5c0", "#fa9fb5", "#f768a1", "#c51b8a", "#7a0177" };
        case 7:
            return { "#feebe2", "#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177" };
        case 8:
            return { "#fff7f3", "#fde0dd", "#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177" };
        case 9:
            return { "#fff7f3", "#fde0dd", "#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177", "#49006a" };
        }
    }
    else if (colorName == "PuRd") {
        switch (colorCount) {
        case 3:
            return { "#e7e1ef", "#c994c7", "#dd1c77" };
        case 4:
            return { "#f1eef6", "#d7b5d8", "#df65b0", "#ce1256" };
        case 5:
            return { "#f1eef6", "#d7b5d8", "#df65b0", "#dd1c77", "#980043" };
        case 6:
            return { "#f1eef6", "#d4b9da", "#c994c7", "#df65b0", "#dd1c77", "#980043" };
        case 7:
            return { "#f1eef6", "#d4b9da", "#c994c7", "#df65b0", "#e7298a", "#ce1256", "#91003f" };
        case 8:
            return { "#f7f4f9", "#e7e1ef", "#d4b9da", "#c994c7", "#df65b0", "#e7298a", "#ce1256", "#91003f" };
        case 9:
            return { "#f7f4f9", "#e7e1ef", "#d4b9da", "#c994c7", "#df65b0", "#e7298a", "#ce1256", "#980043", "#67001f" };
        }
    }
    else if (colorName == "OrRd") {
        switch (colorCount) {
        case 3:
            return { "#fee8c8", "#fdbb84", "#e34a33" };
        case 4:
            return { "#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f" };
        case 5:
            return { "#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000" };
        case 6:
            return { "#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#e34a33", "#b30000" };
        case 7:
            return { "#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#990000" };
        case 8:
            return { "#fff7ec", "#fee8c8", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#990000" };
        case 9:
            return { "#fff7ec", "#fee8c8", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#b30000", "#7f0000" };
        }
    }
    else if (colorName == "YlOrRd") {
        switch (colorCount) {
        case 3:
            return { "#ffeda0", "#feb24c", "#f03b20" };
        case 4:
            return { "#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c" };
        case 5:
            return { "#ffffb2", "#fecc5c", "#fd8d3c", "#f03b20", "#bd0026" };
        case 6:
            return { "#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#f03b20", "#bd0026" };
        case 7:
            return { "#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026" };
        case 8:
            return { "#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026" };
        case 9:
            return { "#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026" };
        }
    }
    else if (colorName == "YlOrBr") {
        switch (colorCount) {
        case 3:
            return { "#fff7bc", "#fec44f", "#d95f0e" };
        case 4:
            return { "#ffffd4", "#fed98e", "#fe9929", "#cc4c02" };
        case 5:
            return { "#ffffd4", "#fed98e", "#fe9929", "#d95f0e", "#993404" };
        case 6:
            return { "#ffffd4", "#fee391", "#fec44f", "#fe9929", "#d95f0e", "#993404" };
        case 7:
            return { "#ffffd4", "#fee391", "#fec44f", "#fe9929", "#ec7014", "#cc4c02", "#8c2d04" };
        case 8:
            return { "#ffffe5", "#fff7bc", "#fee391", "#fec44f", "#fe9929", "#ec7014", "#cc4c02", "#8c2d04" };
        case 9:
            return { "#ffffe5", "#fff7bc", "#fee391", "#fec44f", "#fe9929", "#ec7014", "#cc4c02", "#993404", "#662506" };
        }
    }
    else if (colorName == "Purples") {
        switch (colorCount) {
        case 3:
            return { "#efedf5", "#bcbddc", "#756bb1" };
        case 4:
            return { "#f2f0f7", "#cbc9e2", "#9e9ac8", "#6a51a3" };
        case 5:
            return { "#f2f0f7", "#cbc9e2", "#9e9ac8", "#756bb1", "#54278f" };
        case 6:
            return { "#f2f0f7", "#dadaeb", "#bcbddc", "#9e9ac8", "#756bb1", "#54278f" };
        case 7:
            return { "#f2f0f7", "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#4a1486" };
        case 8:
            return { "#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#4a1486" };
        case 9:
            return { "#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#54278f", "#3f007d" };
        }
    }
    else if (colorName == "Blues") {
        switch (colorCount) {
        case 3:
            return { "#deebf7", "#9ecae1", "#3182bd" };
        case 4:
            return { "#eff3ff", "#bdd7e7", "#6baed6", "#2171b5" };
        case 5:
            return { "#eff3ff", "#bdd7e7", "#6baed6", "#3182bd", "#08519c" };
        case 6:
            return { "#eff3ff", "#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c" };
        case 7:
            return { "#eff3ff", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#084594" };
        case 8:
            return { "#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#084594" };
        case 9:
            return { "#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b" };
        }
    }
    else if (colorName == "Greens") {
        switch (colorCount) {
        case 3:
            return { "#e5f5e0", "#a1d99b", "#31a354" };
        case 4:
            return { "#edf8e9", "#bae4b3", "#74c476", "#238b45" };
        case 5:
            return { "#edf8e9", "#bae4b3", "#74c476", "#31a354", "#006d2c" };
        case 6:
            return { "#edf8e9", "#c7e9c0", "#a1d99b", "#74c476", "#31a354", "#006d2c" };
        case 7:
            return { "#edf8e9", "#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45", "#005a32" };
        case 8:
            return { "#f7fcf5", "#e5f5e0", "#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45", "#005a32" };
        case 9:
            return { "#f7fcf5", "#e5f5e0", "#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45", "#006d2c", "#00441b" };
        }
    }
    else if (colorName == "Oranges") {
        switch (colorCount) {
        case 3:
            return { "#fee6ce", "#fdae6b", "#e6550d" };
        case 4:
            return { "#feedde", "#fdbe85", "#fd8d3c", "#d94701" };
        case 5:
            return { "#feedde", "#fdbe85", "#fd8d3c", "#e6550d", "#a63603" };
        case 6:
            return { "#feedde", "#fdd0a2", "#fdae6b", "#fd8d3c", "#e6550d", "#a63603" };
        case 7:
            return { "#feedde", "#fdd0a2", "#fdae6b", "#fd8d3c", "#f16913", "#d94801", "#8c2d04" };
        case 8:
            return { "#fff5eb", "#fee6ce", "#fdd0a2", "#fdae6b", "#fd8d3c", "#f16913", "#d94801", "#8c2d04" };
        case 9:
            return { "#fff5eb", "#fee6ce", "#fdd0a2", "#fdae6b", "#fd8d3c", "#f16913", "#d94801", "#a63603", "#7f2704" };
        }
    }
    else if (colorName == "Reds") {
        switch (colorCount) {
        case 3:
            return { "#fee0d2", "#fc9272", "#de2d26" };
        case 4:
            return { "#fee5d9", "#fcae91", "#fb6a4a", "#cb181d" };
        case 5:
            return { "#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15" };
        case 6:
            return { "#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15" };
        case 7:
            return { "#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#99000d" };
        case 8:
            return { "#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#99000d" };
        case 9:
            return { "#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15", "#67000d" };
        }
    }
    else if (colorName == "Greys") {
        switch (colorCount) {
        case 3:
            return { "#f0f0f0", "#bdbdbd", "#636363" };
        case 4:
            return { "#f7f7f7", "#cccccc", "#969696", "#525252" };
        case 5:
            return { "#f7f7f7", "#cccccc", "#969696", "#636363", "#252525" };
        case 6:
            return { "#f7f7f7", "#d9d9d9", "#bdbdbd", "#969696", "#636363", "#252525" };
        case 7:
            return { "#f7f7f7", "#d9d9d9", "#bdbdbd", "#969696", "#737373", "#525252", "#252525" };
        case 8:
            return { "#ffffff", "#f0f0f0", "#d9d9d9", "#bdbdbd", "#969696", "#737373", "#525252", "#252525" };
        case 9:
            return { "#ffffff", "#f0f0f0", "#d9d9d9", "#bdbdbd", "#969696", "#737373", "#525252", "#252525", "#000000" };
        }
    }
    else if (colorName == "PuOr") {
        switch (colorCount) {
        case 3:
            return { "#f1a340", "#f7f7f7", "#998ec3" };
        case 4:
            return { "#e66101", "#fdb863", "#b2abd2", "#5e3c99" };
        case 5:
            return { "#e66101", "#fdb863", "#f7f7f7", "#b2abd2", "#5e3c99" };
        case 6:
            return { "#b35806", "#f1a340", "#fee0b6", "#d8daeb", "#998ec3", "#542788" };
        case 7:
            return { "#b35806", "#f1a340", "#fee0b6", "#f7f7f7", "#d8daeb", "#998ec3", "#542788" };
        case 8:
            return { "#b35806", "#e08214", "#fdb863", "#fee0b6", "#d8daeb", "#b2abd2", "#8073ac", "#542788" };
        case 9:
            return { "#b35806", "#e08214", "#fdb863", "#fee0b6", "#f7f7f7", "#d8daeb", "#b2abd2", "#8073ac", "#542788" };
        case 10:
            return { "#7f3b08", "#b35806", "#e08214", "#fdb863", "#fee0b6", "#d8daeb", "#b2abd2", "#8073ac", "#542788", "#2d004b" };
        case 11:
            return { "#7f3b08", "#b35806", "#e08214", "#fdb863", "#fee0b6", "#f7f7f7", "#d8daeb", "#b2abd2", "#8073ac", "#542788", "#2d004b" };
        }
    }
    else if (colorName == "BrBG") {
        switch (colorCount) {
        case 3:
            return { "#d8b365", "#f5f5f5", "#5ab4ac" };
        case 4:
            return { "#a6611a", "#dfc27d", "#80cdc1", "#018571" };
        case 5:
            return { "#a6611a", "#dfc27d", "#f5f5f5", "#80cdc1", "#018571" };
        case 6:
            return { "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e" };
        case 7:
            return { "#8c510a", "#d8b365", "#f6e8c3", "#f5f5f5", "#c7eae5", "#5ab4ac", "#01665e" };
        case 8:
            return { "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#c7eae5", "#80cdc1", "#35978f", "#01665e" };
        case 9:
            return { "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#f5f5f5", "#c7eae5", "#80cdc1", "#35978f", "#01665e" };
        case 10:
            return { "#543005", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#c7eae5", "#80cdc1", "#35978f", "#01665e", "#003c30" };
        case 11:
            return { "#543005", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#f5f5f5", "#c7eae5", "#80cdc1", "#35978f", "#01665e", "#003c30" };
        }
    }
    else if (colorName == "PRGn") {
        switch (colorCount) {
        case 3:
            return { "#af8dc3", "#f7f7f7", "#7fbf7b" };
        case 4:
            return { "#7b3294", "#c2a5cf", "#a6dba0", "#008837" };
        case 5:
            return { "#7b3294", "#c2a5cf", "#f7f7f7", "#a6dba0", "#008837" };
        case 6:
            return { "#762a83", "#af8dc3", "#e7d4e8", "#d9f0d3", "#7fbf7b", "#1b7837" };
        case 7:
            return { "#762a83", "#af8dc3", "#e7d4e8", "#f7f7f7", "#d9f0d3", "#7fbf7b", "#1b7837" };
        case 8:
            return { "#762a83", "#9970ab", "#c2a5cf", "#e7d4e8", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837" };
        case 9:
            return { "#762a83", "#9970ab", "#c2a5cf", "#e7d4e8", "#f7f7f7", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837" };
        case 10:
            return { "#40004b", "#762a83", "#9970ab", "#c2a5cf", "#e7d4e8", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837", "#00441b" };
        case 11:
            return { "#40004b", "#762a83", "#9970ab", "#c2a5cf", "#e7d4e8", "#f7f7f7", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837", "#00441b" };
        }
    }
    else if (colorName == "PiYG") {
        switch (colorCount) {
        case 3:
            return { "#e9a3c9", "#f7f7f7", "#a1d76a" };
        case 4:
            return { "#d01c8b", "#f1b6da", "#b8e186", "#4dac26" };
        case 5:
            return { "#d01c8b", "#f1b6da", "#f7f7f7", "#b8e186", "#4dac26" };
        case 6:
            return { "#c51b7d", "#e9a3c9", "#fde0ef", "#e6f5d0", "#a1d76a", "#4d9221" };
        case 7:
            return { "#c51b7d", "#e9a3c9", "#fde0ef", "#f7f7f7", "#e6f5d0", "#a1d76a", "#4d9221" };
        case 8:
            return { "#c51b7d", "#de77ae", "#f1b6da", "#fde0ef", "#e6f5d0", "#b8e186", "#7fbc41", "#4d9221" };
        case 9:
            return { "#c51b7d", "#de77ae", "#f1b6da", "#fde0ef", "#f7f7f7", "#e6f5d0", "#b8e186", "#7fbc41", "#4d9221" };
        case 10:
            return { "#8e0152", "#c51b7d", "#de77ae", "#f1b6da", "#fde0ef", "#e6f5d0", "#b8e186", "#7fbc41", "#4d9221", "#276419" };
        case 11:
            return { "#8e0152", "#c51b7d", "#de77ae", "#f1b6da", "#fde0ef", "#f7f7f7", "#e6f5d0", "#b8e186", "#7fbc41", "#4d9221", "#276419" };
        }
    }
    else if (colorName == "RdBu") {
        switch (colorCount) {
        case 3:
            return { "#ef8a62", "#f7f7f7", "#67a9cf" };
        case 4:
            return { "#ca0020", "#f4a582", "#92c5de", "#0571b0" };
        case 5:
            return { "#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0" };
        case 6:
            return { "#b2182b", "#ef8a62", "#fddbc7", "#d1e5f0", "#67a9cf", "#2166ac" };
        case 7:
            return { "#b2182b", "#ef8a62", "#fddbc7", "#f7f7f7", "#d1e5f0", "#67a9cf", "#2166ac" };
        case 8:
            return { "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac" };
        case 9:
            return { "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac" };
        case 10:
            return { "#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061" };
        case 11:
            return { "#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061" };
        }
    }
    else if (colorName == "RdGy") {
        switch (colorCount) {
        case 3:
            return { "#ef8a62", "#ffffff", "#999999" };
        case 4:
            return { "#ca0020", "#f4a582", "#bababa", "#404040" };
        case 5:
            return { "#ca0020", "#f4a582", "#ffffff", "#bababa", "#404040" };
        case 6:
            return { "#b2182b", "#ef8a62", "#fddbc7", "#e0e0e0", "#999999", "#4d4d4d" };
        case 7:
            return { "#b2182b", "#ef8a62", "#fddbc7", "#ffffff", "#e0e0e0", "#999999", "#4d4d4d" };
        case 8:
            return { "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#e0e0e0", "#bababa", "#878787", "#4d4d4d" };
        case 9:
            return { "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#ffffff", "#e0e0e0", "#bababa", "#878787", "#4d4d4d" };
        case 10:
            return { "#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#e0e0e0", "#bababa", "#878787", "#4d4d4d", "#1a1a1a" };
        case 11:
            return { "#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#ffffff", "#e0e0e0", "#bababa", "#878787", "#4d4d4d", "#1a1a1a" };
        }
    }
    else if (colorName == "RdYlBu") {
        switch (colorCount) {
        case 3:
            return { "#fc8d59", "#ffffbf", "#91bfdb" };
        case 4:
            return { "#d7191c", "#fdae61", "#abd9e9", "#2c7bb6" };
        case 5:
            return { "#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6" };
        case 6:
            return { "#d73027", "#fc8d59", "#fee090", "#e0f3f8", "#91bfdb", "#4575b4" };
        case 7:
            return { "#d73027", "#fc8d59", "#fee090", "#ffffbf", "#e0f3f8", "#91bfdb", "#4575b4" };
        case 8:
            return { "#d73027", "#f46d43", "#fdae61", "#fee090", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4" };
        case 9:
            return { "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4" };
        case 10:
            return { "#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695" };
        case 11:
            return { "#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695" };
        }
    }
    else if (colorName == "Spectral") {
        switch (colorCount) {
        case 3:
            return { "#fc8d59", "#ffffbf", "#99d594" };
        case 4:
            return { "#d7191c", "#fdae61", "#abdda4", "#2b83ba" };
        case 5:
            return { "#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba" };
        case 6:
            return { "#d53e4f", "#fc8d59", "#fee08b", "#e6f598", "#99d594", "#3288bd" };
        case 7:
            return { "#d53e4f", "#fc8d59", "#fee08b", "#ffffbf", "#e6f598", "#99d594", "#3288bd" };
        case 8:
            return { "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd" };
        case 9:
            return { "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd" };
        case 10:
            return { "#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2" };
        case 11:
            return { "#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2" };
        }
    }
    else if (colorName == "RdYlGn") {
        switch (colorCount) {
        case 3:
            return { "#fc8d59", "#ffffbf", "#91cf60" };
        case 4:
            return { "#d7191c", "#fdae61", "#a6d96a", "#1a9641" };
        case 5:
            return { "#d7191c", "#fdae61", "#ffffbf", "#a6d96a", "#1a9641" };
        case 6:
            return { "#d73027", "#fc8d59", "#fee08b", "#d9ef8b", "#91cf60", "#1a9850" };
        case 7:
            return { "#d73027", "#fc8d59", "#fee08b", "#ffffbf", "#d9ef8b", "#91cf60", "#1a9850" };
        case 8:
            return { "#d73027", "#f46d43", "#fdae61", "#fee08b", "#d9ef8b", "#a6d96a", "#66bd63", "#1a9850" };
        case 9:
            return { "#d73027", "#f46d43", "#fdae61", "#fee08b", "#ffffbf", "#d9ef8b", "#a6d96a", "#66bd63", "#1a9850" };
        case 10:
            return { "#a50026", "#d73027", "#f46d43", "#fdae61", "#fee08b", "#d9ef8b", "#a6d96a", "#66bd63", "#1a9850", "#006837" };
        case 11:
            return { "#a50026", "#d73027", "#f46d43", "#fdae61", "#fee08b", "#ffffbf", "#d9ef8b", "#a6d96a", "#66bd63", "#1a9850", "#006837" };
        }
    }
    else if (colorName == "Accent") {
        switch (colorCount) {
        case 3:
            return { "#7fc97f", "#beaed4", "#fdc086" };
        case 4:
            return { "#7fc97f", "#beaed4", "#fdc086", "#ffff99" };
        case 5:
            return { "#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0" };
        case 6:
            return { "#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#f0027f" };
        case 7:
            return { "#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#f0027f", "#bf5b17" };
        case 8:
            return { "#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#f0027f", "#bf5b17", "#666666" };
        }
    }
    else if (colorName == "Dark2") {
        switch (colorCount) {
        case 3:
            return { "#1b9e77", "#d95f02", "#7570b3" };
        case 4:
            return { "#1b9e77", "#d95f02", "#7570b3", "#e7298a" };
        case 5:
            return { "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e" };
        case 6:
            return { "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02" };
        case 7:
            return { "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d" };
        case 8:
            return { "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666" };
        }
    }
    else if (colorName == "Paired") {
        switch (colorCount) {
        case 3:
            return { "#a6cee3", "#1f78b4", "#b2df8a" };
        case 4:
            return { "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c" };
        case 5:
            return { "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99" };
        case 6:
            return { "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c" };
        case 7:
            return { "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f" };
        case 8:
            return { "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00" };
        case 9:
            return { "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6" };
        case 10:
            return { "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a" };
        case 11:
            return { "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99" };
        case 12:
            return { "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928" };
        }
    }
    else if (colorName == "Pastel1") {
        switch (colorCount) {
        case 3:
            return { "#fbb4ae", "#b3cde3", "#ccebc5" };
        case 4:
            return { "#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4" };
        case 5:
            return { "#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#fed9a6" };
        case 6:
            return { "#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#fed9a6", "#ffffcc" };
        case 7:
            return { "#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#fed9a6", "#ffffcc", "#e5d8bd" };
        case 8:
            return { "#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#fed9a6", "#ffffcc", "#e5d8bd", "#fddaec" };
        case 9:
            return { "#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#fed9a6", "#ffffcc", "#e5d8bd", "#fddaec", "#f2f2f2" };
        }
    }
    else if (colorName == "Pastel2") {
        switch (colorCount) {
        case 3:
            return { "#b3e2cd", "#fdcdac", "#cbd5e8" };
        case 4:
            return { "#b3e2cd", "#fdcdac", "#cbd5e8", "#f4cae4" };
        case 5:
            return { "#b3e2cd", "#fdcdac", "#cbd5e8", "#f4cae4", "#e6f5c9" };
        case 6:
            return { "#b3e2cd", "#fdcdac", "#cbd5e8", "#f4cae4", "#e6f5c9", "#fff2ae" };
        case 7:
            return { "#b3e2cd", "#fdcdac", "#cbd5e8", "#f4cae4", "#e6f5c9", "#fff2ae", "#f1e2cc" };
        case 8:
            return { "#b3e2cd", "#fdcdac", "#cbd5e8", "#f4cae4", "#e6f5c9", "#fff2ae", "#f1e2cc", "#cccccc" };
        }
    }
    else if (colorName == "Set1") {
        switch (colorCount) {
        case 3:
            return { "#e41a1c", "#377eb8", "#4daf4a" };
        case 4:
            return { "#e41a1c", "#377eb8", "#4daf4a", "#984ea3" };
        case 5:
            return { "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00" };
        case 6:
            return { "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33" };
        case 7:
            return { "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628" };
        case 8:
            return { "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf" };
        case 9:
            return { "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999" };
        }
    }
    else if (colorName == "Set2") {
        switch (colorCount) {
        case 3:
            return { "#66c2a5", "#fc8d62", "#8da0cb" };
        case 4:
            return { "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3" };
        case 5:
            return { "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854" };
        case 6:
            return { "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f" };
        case 7:
            return { "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494" };
        case 8:
            return { "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3" };
        }
    }
    else if (colorName == "Set3") {
        switch (colorCount) {
        case 3:
            return { "#8dd3c7", "#ffffb3", "#bebada" };
        case 4:
            return { "#8dd3c7", "#ffffb3", "#bebada", "#fb8072" };
        case 5:
            return { "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3" };
        case 6:
            return { "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462" };
        case 7:
            return { "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69" };
        case 8:
            return { "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5" };
        case 9:
            return { "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9" };
        case 10:
            return { "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd" };
        case 11:
            return { "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5" };
        case 12:
            return { "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f" };
        }
    }
    return {};
}

#endif
