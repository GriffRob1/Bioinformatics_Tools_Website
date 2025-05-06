import React, {useRef, useState} from 'react';
import ToolDescription from './ToolDescription';
import BlueButton from "./BlueButton";


export default function HomeContent({toolsList}) {
    let toolsListRef = useRef(toolsList);
    let tabNumber = useRef(0);

    //assigns toolsListRef to the 9 most popular tools
    toolsListRef.current.sort((tool1,tool2) => {
        if (tool1.popularity > tool2.popularity) return -1;
        if (tool1.popularity < tool2.popularity) return 1;
        else return 0;
    });
    toolsListRef.current = toolsListRef.current.slice(0,10);

    let [toolsListForCurrentTab, setToolsListForCurrentTab] = useState(toolsListRef.current.slice(0,3));

    const tabRightThroughToolsList = () => {
        tabNumber.current = (tabNumber.current + 1) % 3;
        setToolsListForCurrentTab(toolsListRef.current.slice(tabNumber.current * 3, (tabNumber.current * 3) + 3));
    }

    const tabLeftThroughToolsList = () => {
        tabNumber.current = (((tabNumber.current - 1) % 3) + 3) % 3;
        setToolsListForCurrentTab(toolsListRef.current.slice(tabNumber.current * 3, (tabNumber.current * 3) + 3));
    }

    //navigation dots correspond to the tab number
    let navigationDots = [];
    for (let i = 0; i < 3; i++) {
        let dotStyle;
        if (i === tabNumber.current) {
            dotStyle = {backgroundColor: 'black'}
        }
        else {
            dotStyle = {backgroundColor: '#bbb'}
        }
        navigationDots.push(<div key={i} className={'dot'} style={dotStyle}></div>)
    }


    return (
      <section className={'container home-content'}>
          <h2 className={'most-popular-tools-title'}>Most Popular Tools</h2>
          <div className={'container most-popular-tools-container'}>
              <div className={'container triangle-left-container'} onClick={tabLeftThroughToolsList}>
                  <div className={'triangle-left'}></div>
              </div>
              <ul className={'container tools-list'}>
                  {toolsListForCurrentTab.map((tool) => (
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
              <div className={'container triangle-right-container'} onClick={tabRightThroughToolsList}>
                  <div className={'triangle-right'}></div>
              </div>
          </div>
          <div className={'container dot-container'}>{navigationDots}</div>
          <BlueButton content={'See All Tools'}
                      URL={'/tools'}
                      buttonClass={'see-all-tools-button'} />
      </section>
    );
}