import React from 'react';
import {useNavigate} from 'react-router-dom';
import ToolDescription from './ToolDescription';



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


    const navigate = useNavigate();
    const navigateToTools = () => {
        navigate('/tools');
    }

    return (
      <section className={'container home-content'}>
          <h2 className={'most-popular-tools-title'}>Most Popular Tools</h2>
          <div className={'container most-popular-tools-container'}>
              <div className={'triangle-left'}></div>
              <ul className={'container tools-list'}>
                  {toolsSampleData.map((item) => (
                      <li key={item.id}>
                          <ToolDescription URL={item.URL} imagePath={item.imagePath} toolTitle={item.toolTitle} textDescription={item.textDescription}/>
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
          <button className={'blue-button see-all-tools-button'} onClick={navigateToTools}>See All Tools</button>
      </section>
    );
}