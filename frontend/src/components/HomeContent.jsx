import React from 'react';
import ToolDescription from './ToolDescription';
import BlueButton from "./BlueButton";


export default function HomeContent() {

    let toolsSampleData = `[
        {
            "id": 0,
            "popularity": 9,
            "URL": "/genome-assembler",
            "imagePath": "/images/favicon.ico",
            "toolTitle": "Genome Assembler",
            "textDescription": "Assemble a genome using short nucleotide reads"
        },

        {
            "id": 1,
            "popularity": 5,
            "URL": "/motif-finder",
            "imagePath": "/images/logo192.png",
            "toolTitle": "Motif Finder",
            "textDescription": "Find a common pattern in a list of reads"
        },

        {
            "id": 2,
            "popularity": 3,
            "URL": "/sequence-alignment",
            "imagePath": "/images/logo512.png",
            "toolTitle": "Sequence Alignment",
            "textDescription": "Align 2 or more reads using an alignment algorithm"
        }
    ]`

    toolsSampleData = JSON.parse(toolsSampleData)

    return (
      <section className={'container home-content'}>
          <h2 className={'most-popular-tools-title'}>Most Popular Tools</h2>
          <div className={'container most-popular-tools-container'}>
              <div className={'triangle-left'}></div>
              <ul className={'container tools-list'}>
                  {toolsSampleData.map((tool) => (
                      <li key={tool.id}>
                          <ToolDescription
                              id={tool.id}
                              URL={tool.URL}
                              imagePath={tool.imagePath}
                              toolTitle={tool.toolTitle}
                              textDescription={tool.textDescription}
                              popularity={tool.popularity}
                              dateAdded={tool.dateAdded}
                              isFavoritable={false}/>
                      </li>
                  ))}
              </ul>
              <div className={'triangle-right'}></div>
          </div>
          <div className={'container dot-container'}>
              <div className={'dot'}></div>
              <div className={'dot'}></div>
              <div className={'dot'}></div>
          </div>
          <BlueButton content={'See All Tools'}
                      URL={'/tools'}
                      buttonClass={'see-all-tools-button'} />
      </section>
    );
}