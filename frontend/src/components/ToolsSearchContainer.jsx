import React, {useContext} from 'react';
import ToolDescription from "./ToolDescription";
import {AppContext} from "../App";

export default function ToolsSearchContainer({inputName,
                                              inputCategory,
                                              inputSortBy,
                                              inputShowOnlyFavorites}) {

    const toolsList = useContext(AppContext);

    let favoritesList = localStorage.getItem('favorites')
    if (favoritesList) {
        favoritesList = JSON.parse(favoritesList);
    }

    //filters tools by name, category, and favorites
    let filteredTools = [];
    toolsList.forEach((tool) => {
        let inputNameLowerCase = inputName.toLowerCase();
        let toolTitleLowerCase = tool.toolTitle.toLowerCase()
        if (toolTitleLowerCase.indexOf(inputNameLowerCase) === -1) {
            return;
        }
        if (tool.category.indexOf(inputCategory) === -1 ) { //indexOf ensures the default option "" of inputCategory shows all tools
            return;
        }
        if (inputShowOnlyFavorites && (favoritesList.indexOf(tool._id.$oid) === -1)) {
            return;
        }

        filteredTools.push(tool);
    });

    //sorts tools based on user input
    if (inputSortBy === 'sort-by-alphabetical') {
        filteredTools.sort((tool1, tool2) => {
            return (tool1.toolTitle.localeCompare(tool2.toolTitle));
        });
    }
    else if (inputSortBy === 'sort-by-most-popular') {
        filteredTools.sort((tool1,tool2) => {
            if (tool1.popularity > tool2.popularity) return -1;
            if (tool1.popularity < tool2.popularity) return 1;
            else return 0;
        });
    }
    else if (inputSortBy === 'sort-by-newest') {
        filteredTools.sort((tool1,tool2) => {
            const date1 = new Date(tool1.dateAdded)
            const date2 = new Date(tool2.dateAdded)
            if (date1 > date2) return -1;
            if (date1 < date2) return 1;
            else return 0;
        });
    }

    return (
      <section className={'container tools-search-container'}>
          <div className={'grid-container tools-search-container-list'}>
              {filteredTools.map((tool) => (
                  <div key={tool._id.$oid} className={'container tool-description-container'}>
                      <ToolDescription
                          oid={tool._id.$oid}
                          URL={tool.URL}
                          imagePath={tool.imagePath}
                          toolTitle={tool.toolTitle}
                          textDescription={tool.textDescription}
                          popularity={tool.popularity}
                          dateAdded={tool.dateAdded}
                          category={tool.category}/>
                  </div>
              ))}
          </div>
      </section>
    );
}