import React from 'react';
import ToolDescription from "./ToolDescription";

export default function ToolsSearchContainer({toolsSampleData,
                                              inputName,
                                              inputCategory,
                                              inputSortBy}) {

    toolsSampleData = JSON.parse(toolsSampleData)

    let filteredTools = []

    toolsSampleData.forEach((tool) => {
        let inputNameLowerCase = inputName.toLowerCase();
        let toolTitleLowerCase = tool.toolTitle.toLowerCase()
        if (toolTitleLowerCase.indexOf(inputNameLowerCase) === -1) {
            return;
        }

        let inputCategoryLowerCase = inputCategory.toLowerCase();
        let categoryLowerCase = tool.category.toLowerCase();
        if (categoryLowerCase.indexOf(inputCategoryLowerCase) === -1 ) {
            return;
        }

        filteredTools.push(tool);
    });

    if (inputSortBy === 'sort-by-alphabetical') {
        filteredTools.sort((tool1, tool2) => {
            return (tool1.toolTitle.localeCompare(tool2.toolTitle));
        });
    }
    else if (inputSortBy === 'sort-by-most-popular') {
        filteredTools.sort((tool1,tool2) => {
            if (tool1.popularity > tool2.popularity) return -1;
            if (tool1.popularity < tool2.popularity) return 1;
            if (tool1.popularity === tool2.popularity) return 0;
        });
    }
    else if (inputSortBy === 'sort-by-newest') {
        filteredTools.sort((tool1,tool2) => {
            const date1 = new Date(tool1.dateAdded)
            const date2 = new Date(tool2.dateAdded)
            console.log(date1)
            console.log(date2)
            if (date1 < date2) return -1;
            if (date1 > date2) return 1;
            if (date1 === date2) return 0;
        });
    }

    return (
      <section className={'container tools-search-container'}>
          <div className={'grid-container tools-search-container-list'}>
              {filteredTools.map((item) => (
                  <div key={item.id} className={'container tool-description-container'}>
                      <ToolDescription
                          URL={item.URL}
                          imagePath={item.imagePath}
                          toolTitle={item.toolTitle}
                          textDescription={item.textDescription}
                          popularity={item.popularity}
                          dateAdded={item.dateAdded}
                          category={item.category}/>
                  </div>
              ))}
          </div>
      </section>
    );
}